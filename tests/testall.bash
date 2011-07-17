#!/bin/bash
function printtable {
cat >> testall.html <<EOF
	<table>
	<tr><th>Test</hd><th>Error (limit)</hd><th>Runtime [s]</hd><th>Pass</th></tr>
EOF
make test -k -C $1
awk -Ft 'BEGIN{ORS=""}{print "<tr><td>", $7, "</td><td>", $3, " (",$4, ")</td><td>", $6, "</td><td bgcolor=#"; if($5==1){print "00FF00";}else{print "FF0000";} print">", $5 ,"</td></tr>\n";}' $1/error.txt >> testall.html
cat >> testall.html <<EOF
	</table>
EOF

}

cat > testall.html <<EOF
<html>
<head>
	<title>Nbody Unit Testing Summary</title>
</head>
<body>
	<h1>Nbody Unit Testing Summary</h1>
	<p>This file was created automatically by the file 'tests/testall.bash' on `date`. The machine was running `uname -a`. </p>
EOF



cat >> testall.html <<EOF
	<h2>Wisdom Holman Integrator</h2>
EOF
printtable wisdom_holman

cat >> testall.html <<EOF
	<h2>Shearing Sheat (dummy)</h2>
EOF
printtable shearing_sheet





cat >> testall.html <<EOF
</body>
</html>
EOF
