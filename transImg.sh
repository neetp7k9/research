( set -x ; for f_png in *.gif ; do f="${f_png%.gif}" ; convert "$f_png" "$f.pnm" && ../../potrace/potrace "$f.pnm" -s -o "$f.svg" ; done )
