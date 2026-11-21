import sys

with open('docs/index.md', 'r') as f:
    lines = f.readlines()

with open('docs/contributors.txt', 'r') as f:
    snippet = f.read()

new_content = []
skip = False

for line in lines:
    if "<!-- START CONTRIBUTORS -->" in line:
        new_content.append(line)
        new_content.append(snippet + "\n")
        skip = True
    elif "<!-- END CONTRIBUTORS -->" in line:
        new_content.append(line)
        skip = False
    elif not skip:
        new_content.append(line)

with open('docs/index.md', 'w') as f:
    f.writelines(new_content)
