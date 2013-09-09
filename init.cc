/*
 * init.cc - program initialization
 *
 * Copyright © 2013 H.-P. Schadler  <hanspeter.schadler@uni-graz.at>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#include <getopt.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "init.h"
#include "version.h"

int Init(int &argc, char *argv[], GlobalSettings &settings) {
  const char texthelp[]="Usage: lqcd.x [OPTIONS] ... \n"
      "\n"
      "Mandatory arguments to long options are mandatory for short options too.\n"
      "  -s, --Ns SSIZE             spatial lattice extent (default = 4)\n"
      "  -t, --Nt TSIZE             temporal lattice extent (default = 4)\n"
      "  -n, --nmeas NMEAS          number of configurations (default = 1)\n"
      "  -k, --nskip NSKIP          number of skipped configurations (default = 0)\n"
      "  -e, --nequi NEQUI          number of equilibration sweeps (default = 0)\n"
      "  --3d                       write 3dcluster data files\n"
      "  -m, --memory               drop clusterdata vectors after calculating the observables\n"
      "\n"  
      "  -h  --help                 display this help and exit\n"
      "  -v  --version              output version information and exit\n"
      "\n"
      "Exit status:\n"
      " 0  if OK,\n"
      " 1  if minor problems,\n"
      " 2  if serious trouble.\n"
      "\n"
      "Report bugs to hanspeter.schadler@uni-graz.at";

  std::cout << "lqcd.x " << MAJOR_VERSION << "." << MINOR_VERSION << "." 
    << REVISION_VERSION << " ~ " << __DATE__ << " " << __TIME__ 
    << std::endl << std::endl
    << "Lattice QCD Monte Carlo simulation - modern C++ implementation" 
    << std::endl
    << "(C) 2013 Hans-Peter Schadler <hanspeter.schadler@uni-graz.at>"
    << std::endl << std::endl;
  
  if (argc<1) {
    std::cout << texthelp << std::endl;
    return 2;
  }

  std::cout << std::endl << "Initializing... " << std::endl << std::endl;

  // Parse command line options and initialize global settings
  // We first set the default options for lattice size, measurements
  // and other important variables
  
  settings.ns=4;
  settings.nt=4;
  settings.nequi=10;
  settings.nmeas=10;
  settings.nskip=1;
  
  int c;

  while (1){
	  static struct option long_options[] =
		  {
		  /* These options don't set a flag.
		  We distinguish them by their indices. */
		  {"Ns", required_argument, 0, 's'},
		  {"Nt", required_argument, 0, 't'},
		  {"nmeas", required_argument, 0, 'n'},
		  {"nskip", required_argument, 0, 'k'},
		  {"nequi", required_argument, 0, 'e'},
		  {"3d", no_argument, 0, 0},
		  /* These options set a flag. */
		  {"help", no_argument, 0, 'h'},
		  {"version", no_argument, 0, 'v'},
		  {0, 0, 0, 0}
		  };
		
	  /* getopt_long stores the option index here. */
	  int option_index = 0;

	  c = getopt_long (argc, argv, "s:t:n:k:e:hv",
	  long_options, &option_index);

	  /* Detect the end of the options. */
	  if (c == -1) break;

	  switch (c) {
		  case 0:
			  /* If this option set a flag, do nothing else now. */
			  /*if (long_options[option_index].flag != 0)
				  break;
			  printf ("option %s", long_options[option_index].name);
			  if (optarg)
				  printf (" with arg %s", optarg);
			  printf ("\n");A */
	      //            if( strcmp( "3d", long_options[option_index].name ) == 0 )
		    //                 do3d = true;
			  break;

		  case 's':
		    settings.ns = std::atoi(optarg);
			  break;

		  case 't':
		    settings.nt = std::atoi(optarg);
			  break;
			  
		  case 'n':
			  settings.nmeas = std::atoi(optarg);
			  break;

		  case 'k':
			  settings.nskip = std::atoi(optarg);
			  break;
		  
      case 'e':
			  settings.nequi = std::atoi(optarg);
			  break;

		  case 'v':
        std::cout << std::endl 
          << "lqcd.x " << MAJOR_VERSION << "." << MINOR_VERSION << "." 
          << REVISION_VERSION << " ~ " << __DATE__ << " " << __TIME__ 
          << std::endl << std::endl
          << "Lattice QCD Monte Carlo simulation - modern C++ implementation" 
          << std::endl
          << "(C) 2013 Hans-Peter Schadler <hanspeter.schadler@uni-graz.at>"
          << std::endl << std::endl;
			  exit(0);

		  case 'h':
			  std::cout << texthelp << std::endl;
			  exit(0);

		  default:
			  std::cout << texthelp << std::endl;
			  exit(0);
	  }
  }
  
  /*
  string finname;
  // Print any remaining command line arguments (not options).
  if (optind < argc)
  {
	  finname = argv[optind];
  }else{
	  cout << "ERROR: No configuration file specified!" << endl;
  }
  */

  // Perform further initializations
  settings.dim=4;
  settings.nsites=settings.ns*settings.ns*settings.ns*settings.nt;
  
  return 0;
}
