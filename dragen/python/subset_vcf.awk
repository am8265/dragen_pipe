#!/usr/bin/awk -f


BEGIN{
    FS="\t"
    i=0
}

{
    if($0 ~/\#/) print; else if ($7 == "PASS") print
}

END{}
