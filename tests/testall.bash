#!/bin/bash
function printtable {
cat >> testall.html <<EOF
	<table width="100%">
	<tr><th align=left>Test</hd><th align=left>Error (limit)</hd><th align=left>Runtime [s]</hd><th align=left>Pass</th></tr>
EOF
make test -k -C $1
awk -Ft 'BEGIN{ORS="";FS="\t"}{print "<tr><td>", $7, "</td><td width="300px">", $3, " (",$4, ")</td><td width="150px">", $6, "</td><td width="100px" bgcolor=#"; if($5==1){print "00FF00";}else{print "FF0000";} print">", $5 ,"</td></tr>\n";}' $1/error.txt >> testall.html
cat >> testall.html <<EOF
	</table>
EOF
}

function printimage {
cat >> testall.html <<EOF
	<p><img src="$1" /></p>
EOF
}

rm ./*/error.txt

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
	<h2>Accuracy integrator wh</h2>
EOF
printtable accuracy_integrator_wh

cat >> testall.html <<EOF
	<h2>Accuracy tree force</h2>
EOF
printtable accuracy_tree_force
printimage accuracy_tree_force/plot1.png

cat >> testall.html <<EOF
	<h2>Accuracy tree energy</h2>
EOF
printtable accuracy_tree_energy

cat >> testall.html <<EOF
	<h2>Scaling tree</h2>
EOF
make test -k -C scaling_tree
printimage scaling_tree/plot.png


cat >> testall.html <<EOF
	<h2>Shearing Sheat (gravity tree, collision tree)</h2>
EOF
printtable speed_shearing_sheet

cat >> testall.html <<EOF
	<h2>Shearing Sheat (dummy)</h2>
EOF
printtable shearing_sheet

cat >> testall.html <<EOF
	<h2>Bouncing String (collision detection with tree)</h2>
EOF
printtable bouncing_string_tree





cat >> testall.html <<EOF
</body>
</html>
EOF
