#!/bin/bash
if [ -f documentation_lock.txt ]; then
    echo "Documentation update already in progress. Exiting."
    exit
fi
touch documentation_lock.txt

echo "Updating REBOUND"
echo "----------------"
git stash
git pull
repository_head=` git rev-parse HEAD`
documentation_head=$(<documentation_head.txt)
echo "Repository head:    $repository_head"
echo "Documentation head: $documentation_head"

if [[ "$repository_head" == "$documentation_head" ]]; then
    echo "Documentation is up-to-date."
    rm documentation_lock.txt
    git stash pop
    exit
fi

echo "Installing REBOUND"
echo "------------------"
pip install -e .

echo "Generating integrator docs"
echo "--------------------------"
python docs/make_integrators.py

echo "Generating C examples"
echo "---------------------"
python docs/make_c_examples.py

echo "Running mkdocs"
echo "--------------"
mkdocs build  

echo "Documentation updated. Releasing lock"
echo "-------------------------------------"
echo "$repository_head" > documentation_head.txt
rm documentation_lock.txt
git reset --hard HEAD
git stash pop

exit
