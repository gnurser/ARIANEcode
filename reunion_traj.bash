#!/bin/bash
# options with colons afterwards require argiuments; those without don't. See
# http://wiki.bash-hackers.org/howto/getopts_tutorial
# http://tuxtweaks.com/2014/05/bash-getopts/
orca_grid=
multiyear=False
public=False
data_root=~/DATA/NEMO/ARIANE/MH370/Reunion_MH370_mob5/ARIANE_OUTPUT
while getopts r:o:pm FLAG; do
    case $FLAG in
	r)  #get run name
	    runname=$OPTARG
	    #echo "-a used: $OPTARG"
	    #echo "package = $package"
	    ;;
	p) #unrecognized option
	    public=True
	    ;;
	o) #unrecognized option
	    orca_grid=$OPTARG
	    ;;
	m) #unrecognized option
	    multiyear=True
	    ;;
	\?) #unrecognized option
	    exit 2
	    ;;
    esac
done

#This tells getopts to move on to the next argument, so now $1 etc
# relate to further arguments noy hadled by getopts
shift $((OPTIND-1))
### end of getopts code ###

echo "runname = $runname"
echo "orca_grid = $orca_grid"

if [ "$orca_grid" == "4" ]; then
    data_dir=${data_root}/025
elif [ "$orca_grid" == "12" ]; then
    data_dir=${data_root}/1_12
else
    data_dir=${data_root}/1_12_forwards
fi

if [ "$multiyear" != "False" ]; then
    for year in 2006 2007 2008 2009 2010; do
    #for year in 2006 2007; do
	echo "doing year $year"
	dirname=${runname}_${orca_grid}_$year
	yearm1=$(( $year - 1 ))
	if [ ! -d $dirname ]; then
	    mkdir $dirname
	fi
	cd $dirname
	echo "In directory $PWD"

	if [ "$runname" == "fewdots" ]; then
	    traj_movie.py --movie --ntmax 147 -r $data_dir --dotsize 12 --dotalpha 0.5 \
			  0pc_$year 1pc_$year 2pc_$year --date --end_date $year-07-29 --crash_date $yearm1-03-08
	elif [ "$runname" = "lotsdots" ]; then
	    traj_movie.py --movie --ntmax 102 -r $data_dir --dotsize 3 --dotalpha 0.3 \
			  0pc_lots$year 1pc_lots$year 2pc_lots$year 3pc_lots$year \
			  --date --end_date $year-07-29 --crash_date $yearm1-03-08
	elif [ "$runname" = "lotsdotsS" ]; then
	    traj_movie.py --movie --ntmax 102 -r $data_dir --dotsize 3 --dotalpha 0.3 \
			  0pc_lots$year 1pc_lots$year 2pc_lots$year 3pc_lots$year \
			  --date --end_date $year-07-29 --crash_date $yearm1-03-08 \
			  --subarea 80 -36 100 -36 100 -39 80 -39
	elif [ "$runname" = "lotsdotsN" ]; then
	    traj_movie.py --movie --ntmax 102 -r $data_dir --dotsize 3 --dotalpha 0.3 \
			  0pc_lots$year 1pc_lots$year 2pc_lots$year 3pc_lots$year \
			  --date --end_date $year-07-29 --crash_date $yearm1-03-08 \
			  --subarea 80 -36 100 -36 100 -32 80 -32
	    # traj_movie.py --movie --ntmax 102 -r $data_dir --dotsize 3 --dotalpha 0.3 2pc_lots$year --date --end_date $year-07-29 --crash_date $yearm1-03-08
	fi
	cd ..
    done
fi

if [ "$public" == "True" ]; then
    for year in 2010; do
	yearm1=$(( $year - 1 ))
	if [ ! -d public ]; then
	    mkdir public
	fi
	cd public
	if [ "$runname" == "fewdots" ]; then
	    traj_movie.py --movie --ntmax 102 -r $data_dir --dotsize 12 --dotalpha 0.5 -0pc_$year 1pc_$year 2pc_$year --end_date $year-07-29 --crash_date $yearm1-03-08 --movie
	elif [ "$runname" = "lotsdots" ]; then
	    traj_movie.py --movie --ntmax 102 -r $data_dir --dotsize 3 --dotalpha 0.3 0pc_lots$year 1pc_lots$year 2pc_lots$year 3pc_lots$year --end_date $year-07-29 --crash_date $yearm1-03-08
	fi
	cd ..
    done
fi
