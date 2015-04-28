#include <stdio.h>
#include <stdlib.h>

#define BUFFER_SIZE 1000000

char usage[] = "a.out xxx.fa" ;

char buffer[ BUFFER_SIZE + 1 ] ;

int main( int argc, char *argv[] )
{
	FILE *fp = fopen( argv[1], "r" ) ;
	int i, j, k ;
	int bs, l ;
	int state ; // 0-header state, 1-sequence
	
	if ( argc < 2 )
	{
		printf( "%s\n", usage ) ;
		exit( 1 ) ;
	}
	bs = -1 ;
	i = 0 ;
	l = 0 ;
	state = 0 ;
	while ( 1 )	
	{
		if ( i >= bs )
		{
			if ( ( bs = fread( buffer, sizeof( char ), BUFFER_SIZE, fp ) ) == 0 )
				break ;

			i = 0 ;
		}
		if ( buffer[i] == '>' )
		{
			if ( l > 0 )
				printf( "\n" ) ;

			state = 0 ;
			l = 0 ;
		}
		else if ( buffer[i] == '\n' && state == 0 )
		{
			state = 1 ;
			l = -1 ;
		}
		else if ( state == 1 && ( buffer[i] == 'N' || buffer[i] == 'n' || buffer[i] == '\n' ) )
		{
			//printf( "hi\n" ) ;
			++i ;
			continue ;
		}
		
		printf( "%c", buffer[i] ) ;
		++l ;
		if ( l == 80 && state == 1 )
		{
			printf( "\n" ) ;
			l = 0 ;
		}

		++i ;
	}
	return 0 ;
}
