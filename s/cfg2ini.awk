{gsub(/ *$/,"",$3); gsub(/ *$/,"",$1); if ($1!="XTERM") {  print "LF_" $1 "=\"" $3 "\""}}
