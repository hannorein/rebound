#!/bin/bash
# This script provides some consistent linting of all .c and .h files.

# All other files
for filename in src/*.c src/*.h; do
    echo $filename
    ex -u NONE -c ":filetype plugin indent on |:set cinoptions+=+0 |:set smartindent |:set tabstop=4 |:set shiftwidth=4 |:set expandtab | normal! gg=G" -cwq $filename
done

