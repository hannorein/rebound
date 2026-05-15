#!/bin/bash

git log --format='%aN' | sort -u | while read name; do \
    echo -en "$name,"; \
    git log --author="$name" --pretty=tformat: --numstat | \
    awk '{ add += $1; subs += $2 } END { printf "%s,", add }'; \
    git rev-list --count --author="$name" HEAD        
done > contributors.txt

# Manually clean up author aliases
sed \
    -e 's/shangfei/Shangfei Liu/g' \
    -e 's/Liu Shangfei/Shangfei Liu/g' \
    -e 's/pete bartram/Peter Bartram/g' \
    -e 's/Pete Bartram/Peter Bartram/g' \
    -e 's/PeterBartram/Peter Bartram/g' \
    -e 's/Ruth-Huang6012/Ruth Huang/g' \
    -e 's/Ruth-Huang/Ruth Huang/g' \
    -e 's/Hanno REIN/Hanno Rein/g' \
    -e 's/Dan Tamayo/Daniel Tamayo/g' \
    -e 's/dtamayo/Daniel Tamayo/g' \
    -e 's/silburt/Ari Silburt/g' \
    -e 's/Garett\,/Garett Brown\,/g' \
    -e 's/Garett Brown Brown/Garett Brown/g' \
    -e 's/rmelikyan/Robert Melikyan/g' \
    -e 's/tigerchenlu98/Tiger Lu/g' \
    -e 's/Dave\,/David Spiegel\,/g' \
    contributors.txt > contributors_clean.txt

echo "| Name | Number of contributed lines | Number of contributed commits |" >  docs/contributors.txt
echo "| ---- | --------------------------- | ----------------------------- |" >> docs/contributors.txt

awk -F',' '
BEGIN { OFS="," } { added[$1] += $2; commits[$1] += $3 }
END {
    for (name in added) {
        printf "| %s | %s | %s |\n", name,  added[name], commits[name]
    }
}' contributors_clean.txt | sort -rnk3 -t\| | head -n 20 >> docs/contributors.txt
rm contributors_clean.txt
rm contributors.txt

