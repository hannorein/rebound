#!/bin/bash
for filename in src/*.c src/*.h; do
    echo $filename
    ex -u NONE -c ":filetype plugin indent on |:set smartindent |:set tabstop=4 |:set shiftwidth=4 |:set expandtab | normal! gg=G" -cwq $filename
done

