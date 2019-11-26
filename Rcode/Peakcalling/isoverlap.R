# This script used to compare max(A.start,B.start)<=min(A.end,B.end),which means overlap or not
args=commandArgs(T)
max(as.numeric(args[1]),as.numeric(args[2]))<=min(as.numeric(args[3]),as.numeric(args[4])) 
