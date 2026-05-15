#!/bin/bash

git log --format='%aN' | sort -u | while read name; do \
    echo -en "$name,"; \
    git log --author="$name" --pretty=tformat: --numstat | \
    awk '{ add += $1; subs += $2 } END { printf "%s,", add }'; \
    git rev-list --count --author="$name" HEAD        
done > contributors.txt

sed \
    -e 's/shangfei/Shangfei Liu/g' \
    -e 's/Liu Shangfei/Shangfei Liu/g' \
    -e 's/pete bartram/Peter Bartram/g' \
    -e 's/Pete Bartram/Peter Bartram/g' \
    -e 's/PeterBartram/Peter Bartram/g' \
    -e 's/Ruth-Huang/Ruth Huag/g' \
    -e 's/Hanno REIN/Hanno Rein/g' \
    -e 's/Dan Tamayo/Daniel Tamayo/g' \
    -e 's/dtamayo/Daniel Tamayo/g' \
    -e 's/silburt/Ari Silburt/g' \
    -e 's/Garett/Garett Brown/g' \
    contributors.txt > contributors_clean.txt

awk -F',' '
BEGIN { OFS="," } { added[$1] += $2; commits[$1] += $3 }
END {
    for (name in added) {
        print name, added[name], commits[name]
    }
}' contributors_clean.txt | sort -rnk3 -t, > contributors_clean2.txt

#| sort -rnk5 | head -n 10
