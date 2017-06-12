#!/bin/bash -l

body='{
"request": {
"branch":"master"
}}'

curl -s -X POST \
   -H "Content-Type: application/json" \
   -H "Accept: application/json" \
   -H "Travis-API-Version: 3" \
   -H "Authorization: token "$TRAVIS_API_TOKEN \
   -d "$body" \
   https://api.travis-ci.org/repo/dtamayo%2Freboundx/requests
