#!/bin/bash

function check_or_mkdir {
    echo -n "Creating $1 ... "
    if [[ -d $1 && ! -n `find $1 -prune -empty -type d` ]]; then
        echo "Directory exists - skipping it!"
        return `false`
    else 
        echo "Done"
        mkdir -p $1
        return `true`
    fi
}

function check_or_mkdir_no_fail {
    echo -n "Creating $1 ... "
    if [[ -d $1 && ! -n `find $1 -prune -empty -type d` ]]; then
        echo "Directory exists already! Continuing"
        return `true`
    else 
        echo "Done"
        mkdir -p $1
        return `true`
    fi
}



## Functions
function validate_url(){
  if [[ `wget --reject="index.html*" -S --spider $1  2>&1 | egrep 'HTTP/1.1 200 OK|File .* exists.'` ]]; then echo "true"; fi
}
export -f validate_url

function c_echo() {
        printf "\033[34m$*\033[0m\n"
}

progressfilt () {
    # from http://stackoverflow.com/a/4687912/299878
    local flag=false c count cr=$'\r' nl=$'\n'
    while IFS='' read -d '' -rn 1 c
    do
        if $flag
        then
            printf '%c' "$c"
        else
            if [[ $c != $cr && $c != $nl ]]
            then
                count=0
            else
                ((count++))
                if ((count > 1))
                then
                    flag=true
                fi
            fi
        fi
    done
}
