# function defined
fnDOY2yyyymmdd()
{
    date="$1"    # assume all dates are in YYYYMMM format
     echo $date
    year=${date%???}
     echo $year
    jday=${date#$year}
    echo $jday
    for m in `seq 1 12`; do
        for d in `seq 1 31`; do
            yyyymmdd=$(printf "%04d%02d%02d" $year $m $d)
#           echo $yyyymmdd
# 	   doy=$(date +"%j" -d "$yyyymmdd")
     	    #doy=$(date +"%j" -d "$yyyymmdd")
 #             echo $doy
            if test "$jday" = "$doy"; then
                echo "$yyyymmdd"
                return $yyyymmdd
            fi
        done
    done
    echo "Invalid date" >&2
    return 0
}


#calling hello world function
fnDOY2yyyymmdd $1
