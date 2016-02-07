#!/usr/bin/env python

import sys, os, subprocess, inspect
import platform, multiprocessing
import string, re
from datetime import datetime, date, time
import copy
from argparse import ArgumentParser, FileType


def compare_scm(centrifuge_out, true_out):
    db_dic = {}
    first = True
    for line in open(centrifuge_out):
        if first:
            first = False
            continue
        read_name, seq_id, tax_id, score, _, _, _ = line.strip().split('\t')
        if read_name not in db_dic:
            db_dic[read_name] = []
        db_dic[read_name].append(tax_id)

    classified, unclassified, unique_classified = 0, 0, 0
    for line in open(true_out):
        if line.startswith('@'):
            continue
        
        read_name, tax_id = line.strip().split('\t')[:2]
        if read_name not in db_dic:
            unclassified += 1
            continue

        maps = db_dic[read_name]
        if tax_id in maps:
            classified += 1
            if len(maps) == 1:
                unique_classified += 1
        else:
            unclassified += 1

    raw_unique_classified = 0
    for value in db_dic.values():
        if len(value) == 1:
            raw_unique_classified += 1
    return classified, unique_classified, unclassified, len(db_dic), raw_unique_classified


def compare_abundance(centrifuge_out, true_out, debug):
    db_dic = {}
    first = True
    for line in open(centrifuge_out):
        if first:
            first = False
            continue
        genome_name, tax_id, tax_rank, genome_len, num_reads, num_unique_reads, abundance = line.strip().split('\t')
        db_dic[tax_id] = float(abundance)

    SSR = 0.0 # Sum of squared residuals
    first = True
    for line in open(true_out):
        if first:
            first = False
            continue
        
        tax_id, genome_len, num_reads, abundance, genome_name = line.strip().split('\t')
        abundance = float(abundance)
        if tax_id in db_dic:
            SSR += (abundance - db_dic[tax_id]) ** 2;
            if debug:
                print >> sys.stderr, "\t\t\t{:<10}: {:.6} vs. {:.6} (truth vs. centrifuge)".format(tax_id, abundance, db_dic[tax_id])
        else:
            SSR += (abundance) ** 2

    return SSR


def sql_execute(sql_db, sql_query):
    sql_cmd = [
        "sqlite3", sql_db,
        "-separator", "\t",
        "%s;" % sql_query
        ]
    # print >> sys.stderr, sql_cmd
    sql_process = subprocess.Popen(sql_cmd, stdout=subprocess.PIPE)
    output = sql_process.communicate()[0][:-1]
    return output


def create_sql_db(sql_db):
    if os.path.exists(sql_db):
        print >> sys.stderr, sql_db, "already exists!"
        return
    
    columns = [
        ["id", "integer primary key autoincrement"],
        ["index", "text"],
        ["head", "text"],
        ["endType", "text"],
        ["type", "text"],
        ["aligner", "text"],
        ["version", "text"],
        ["numReads", "integer"],
        ["mappedReads", "integer"],
        ["uniqueMappedReads", "integer"],
        ["unmappedReads", "integer"],
        ["time", "real"],
        ["host", "text"],
        ["created", "text"],
        ["cmd", "text"]
        ]
    
    sql_create_table = "CREATE TABLE Classification ("
    for i in range(len(columns)):
        name, type = columns[i]
        if i != 0:
            sql_create_table += ", "
        sql_create_table += ("%s %s" % (name, type))
    sql_create_table += ");"
    sql_execute(sql_db, sql_create_table)


def write_analysis_data(sql_db, genome_name, database_name):
    if not os.path.exists(sql_db):
        return
    
    aligners = []
    sql_aligners = "SELECT aligner FROM ReadCosts GROUP BY aligner"
    output = sql_execute(sql_db, sql_aligners)
    aligners = output.split()

    can_read_types = ["all", "M", "2M_gt_15", "2M_8_15", "2M_1_7", "gt_2M"]    
    tmp_read_types = []
    sql_types = "SELECT type FROM ReadCosts GROUP BY type"
    output = sql_execute(sql_db, sql_types)
    tmp_read_types = output.split()

    read_types = []
    for read_type in can_read_types:
        if read_type in tmp_read_types:
            read_types.append(read_type)

    for paired in [False, True]:
        database_fname = genome_name + "_" + database_name
        if paired:
            end_type = "paired"
            database_fname += "_paired"
        else:
            end_type = "single"
            database_fname += "_single"
        database_fname += ".analysis"
        database_file = open(database_fname, "w")
        print >> database_file, "end_type\ttype\taligner\tnum_reads\ttime\tmapped_reads\tunique_mapped_reads\tunmapped_reads\tmapping_point\ttrue_gtf_junctions\ttemp_junctions\ttemp_gtf_junctions"
        for aligner in aligners:
            for read_type in read_types:
                sql_row = "SELECT end_type, type, aligner, num_reads, time, mapped_reads, unique_mapped_reads, unmapped_reads, mapping_point, true_gtf_junctions, temp_junctions, temp_gtf_junctions FROM ReadCosts"
                sql_row += " WHERE genome = '%s' and head = '%s' and aligner = '%s' and type = '%s' and end_type = '%s' ORDER BY created DESC LIMIT 1" % (genome_name, database_name, aligner, read_type, end_type)
                output = sql_execute(sql_db, sql_row)
                if output:
                    print >> database_file, output

        database_file.close()


def evaluate(index_base,
             num_fragment,
             paired,
             error_rate,
             programs,
             runtime_only,
             sql,
             verbose,
             debug):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(evaluate))
    path_base = os.path.dirname(curr_script)

    sql_db_name = "analysis.db"
    if not os.path.exists(sql_db_name):
        create_sql_db(sql_db_name)

    num_cpus = multiprocessing.cpu_count()
    if num_cpus > 8:
        num_threads = min(8, num_cpus)
        desktop = False
    else:
        num_threads = min(3, num_cpus)
        desktop = True

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    # Check if indexes exists, otherwise create indexes
    index_path = "%s/indexes/Centrifuge" % path_base
    if not os.path.exists(path_base + "/indexes"):
        os.mkdir(path_base + "/indexes")
    if not os.path.exists(index_path):
        os.mkdir(index_path)
    index_fnames = ["%s/%s.%d.cf" % (index_path, index_base, i+1) for i in range(3)]
    assert check_files(index_fnames)
    
        
    # Check if simulated reads exist, otherwise simulate reads
    read_path = "%s/reads" % path_base
    if not os.path.exists(read_path):
        os.mkdir(read_path)
    read_base = "%s_%dM" % (index_base, num_fragment / 1000000)
    if error_rate > 0.0:
        read_base += "%.2fe" % error_rate

    read1_fname = "%s/%s_1.fa" % (read_path, read_base)
    read2_fname = "%s/%s_2.fa" % (read_path, read_base)
    truth_fname = "%s/%s.truth" % (read_path, read_base)
    scm_fname = "%s/%s.scm" % (read_path, read_base)
    read_fnames = [read1_fname, read2_fname, truth_fname, scm_fname]
    if not check_files(read_fnames):
        print >> sys.stderr, "Simulating reads %s_1.fq %s_2.fq ..." % (read_base, read_base)
        centrifuge_simulate = os.path.join(path_base, "centrifuge_simulate_reads.py")
        simulate_cmd = [centrifuge_simulate,
                        "--num-fragment", str(num_fragment)]
        if error_rate > 0.0:
            simulate_cmd += ["--error-rate", str(error_rate)]
        simulate_cmd += ["%s/%s" % (index_path, index_base),
                         "%s/%s" % (read_path, read_base)]
        
        simulate_proc = subprocess.Popen(simulate_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        simulate_proc.communicate()
        assert check_files(read_fnames)

    if runtime_only:
        verbose = True

    if paired:
        base_fname = read_base + "_paired"
    else:
        base_fname = read_base + "_single"

    print >> sys.stderr, "Database: %s" % (index_base)
    if paired:
        print >> sys.stderr, "\t%d million pairs" % (num_fragment / 1000000)
    else:
        print >> sys.stderr, "\t%d million reads" % (num_fragment / 1000000)

    aligner_bin_base = "%s/.." % path_base
    def get_aligner_version(aligner, version):
        version = ""
        if aligner == "centrifuge":
            if version:
                cmd = ["%s/%s_%s/%s" % (aligner_bin_base, aligner, version, aligner)]
            else:
                cmd = ["%s/%s" % (aligner_bin_base, aligner)]
            cmd += ["--version"]                    
            cmd_process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            version = cmd_process.communicate()[0][:-1].split("\n")[0]
            version = version.split()[-1]
        else:
            assert False

        return version

    def get_aligner_cmd(aligner, version, read1_fname, read2_fname, out_fname):
        cmd = []
        if aligner == "centrifuge":
            if version:
                cmd = ["%s/centrifuge_%s/centrifuge" % (aligner_bin_base, version)]
            else:
                cmd = ["%s/centrifuge" % (aligner_bin_base)]
            cmd += ["-f",
                    "-p", str(num_threads),
                    "%s/%s" % (index_path, index_base)]
            # cmd += ["-k", "5"]
            if paired:
                cmd += ["-1", read1_fname,
                        "-2", read2_fname]
            else:
                cmd += ["-U", read1_fname]                        
        else:
            assert False

        return cmd

    init_time = {"centrifuge" : 0.0}
    for program, version in programs:
        program_name = program
        if version:
            program_name += ("_%s" % version)

        print >> sys.stderr, "\t%s\t%s" % (program_name, str(datetime.now()))
        if paired:
            program_dir = program_name + "_paired"
        else:
            program_dir = program_name + "_single"
            
        if not os.path.exists(program_dir):
            os.mkdir(program_dir)
        os.chdir(program_dir)

        out_fname = "centrifuge.output"
        if runtime_only:
            out_fname = "/dev/null"

        duration = -1.0
        if not os.path.exists(out_fname):
            # Classify all reads
            program_cmd = get_aligner_cmd(program, version, read1_fname, read2_fname, out_fname)
            start_time = datetime.now()
            if verbose:
                print >> sys.stderr, "\t", start_time, " ".join(program_cmd)
            if program in ["centrifuge"]:
                proc = subprocess.Popen(program_cmd, stdout=open(out_fname, "w"), stderr=subprocess.PIPE)
            else:
                proc = subprocess.Popen(program_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc.communicate()
            finish_time = datetime.now()
            duration = finish_time - start_time
            assert program in init_time
            duration = duration.total_seconds() - init_time[program]
            if duration < 0.1:
                duration = 0.1
            if verbose:
                print >> sys.stderr, "\t", finish_time, "finished:", duration

        if runtime_only:
            os.chdir("..")
            continue

        classified, unique_classified, unclassified, raw_classified, raw_unique_classified = compare_scm(out_fname, scm_fname)
        assert classified + unclassified == num_fragment

        print >> sys.stderr, "\t\tclassified: {:,}, multi classified: {:,}".format(raw_unique_classified, raw_classified)
        print >> sys.stderr, "\t\tcorrectly classified: {:,} ({:.2}%%)".format(classified, float(classified) * 100.0 / num_fragment)
        print >> sys.stderr, "\t\tuniquely and correctly classified: {:,} ({:.2}%%)".format(unique_classified, float(unique_classified) * 100.0 / num_fragment)

        # Calculate sum of squared residuals in abundance
        abundance_SSR = compare_abundance("centrifuge_report.csv", truth_fname, debug)
        print >> sys.stderr, "\t\tsum of squared residuals in abundance: {}".format(abundance_SSR)

        """
        if duration > 0.0:
            if sql and os.path.exists("../" + sql_db_name):
                if paired:
                    end_type = "paired"
                else:
                    end_type = "single"
                sql_insert = "INSERT INTO \"ReadCosts\" VALUES(NULL, '%s', '%s', '%s', '%s', '%s', '%s', %d, %d, %d, %d, %f, %f, %d, %d, %d, '%s', datetime('now', 'localtime'), '%s');" % \
                    (genome, data_base, end_type, readtype2, aligner_name, get_aligner_version(aligner, version), numreads, mapped, unique_mapped, unmapped, mapping_point, duration, len(junctions), temp_junctions, temp_gtf_junctions, platform.node(), " ".join(aligner_cmd))
                sql_execute("../" + sql_db_name, sql_insert)     

            if two_step:
                align_stat.append([readtype2, aligner_name, numreads, duration, mapped, unique_mapped, unmapped, mapping_point, len(junctions), temp_junctions, temp_gtf_junctions])
            else:
                align_stat[-1].extend([numreads, duration, mapped, unique_mapped, unmapped, mapping_point, len(junctions), temp_junctions, temp_gtf_junctions])
 
        os.system("touch %s.done" % type_sam_fname2)
        """

        os.chdir("..")

        """
        if os.path.exists(sql_db_name):
            write_analysis_data(sql_db_name, genome, data_base)
        """
        

if __name__ == "__main__":
    parser = ArgumentParser(
        description='Centrifuge evaluation')
    parser.add_argument('index_base',
                        nargs='?',
                        type=str,
                        help='Centrifuge index')
    parser.add_argument('--num-fragment',
                        dest="num_fragment",
                        action='store',
                        type=int,
                        default=1,
                        help='Number of fragments in millions (default: 1)')
    parser.add_argument('--paired',
                        dest='paired',
                        action='store_true',
                        help='Paired-end reads')
    parser.add_argument('--error-rate',
                        dest='error_rate',
                        action='store',
                        type=float,
                        default=0.0,
                        help='per-base sequencing error rate (%%) (default: 0.0)')    
    parser.add_argument("--program-list",
                        dest="programs",
                        type=str,
                        default="centrifuge",
                        help="A comma-separated list of aligners (default: centrifuge)")
    parser.add_argument('--runtime-only',
                        dest='runtime_only',
                        action='store_true',
                        help='Just check runtime without evaluation')    
    parser.add_argument('--no-sql',
                        dest='sql',
                        action='store_false',
                        help='Do not write results into a sqlite database')
    parser.add_argument("--simulate-interval",
                        dest="simulate_interval",
                        type=int,
                        default=1,
                        help="Reads simulated at every these base pairs (default: 1)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')
    parser.add_argument('--debug',
                        dest='debug',
                        action='store_true',
                        help='Debug')

    args = parser.parse_args()
    if not args.index_base:
        parser.print_help()
        exit(1)

    programs = []
    for program in args.programs.split(','):
        if '_' in program:
            programs.append(program.split('_'))
        else:
            programs.append([program, ""])
            
    evaluate(args.index_base,
             args.num_fragment * 1000000,
             args.paired,
             args.error_rate,
             programs,
             args.runtime_only,
             args.sql,
             args.verbose,
             args.debug)
