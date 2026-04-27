#!/bin/bash
# This script provides some consistent linting of all .c and .h files.

# Special options for rebound.h
ex -u NONE -c ":filetype plugin indent on |:set cinoptions+=+0 |:set smartindent |:set tabstop=4 |:set shiftwidth=4 |:set expandtab | normal! gg=G" -cwq src/rebound.h

# All other files
for filename in src/*.c src/*.h; do
    if [ "$filename" == "src/rebound.h" ]; then
        continue
    fi
    echo $filename
    ex -u NONE -c ":filetype plugin indent on |:set smartindent |:set tabstop=4 |:set shiftwidth=4 |:set expandtab | normal! gg=G" -cwq $filename
done

