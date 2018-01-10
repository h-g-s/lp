\ENCODING=ISO-8859-1
\Problem name: mip

Minimize
 obj: - v1 + v2
Subject To
 rest1: v1 + v2 >= 3.5
 rest2: - v1 + 5 v2 <= 22.5
 rest3: 2.5 v1 + v2 <= 19.5
 rest4: v1 - 5.5 v2 <= 2.5
 rest5: - v1 + v2 <= 2.5
 rest6: 8.5 v1 + v2 <= 58.5
 rest7: 2.5 v1 + 5.5 v2 >= 12
 rest8: - 1.5 v1 + 5 v2 <= 17.75
Bounds
      v1 >= 0
 2 <= v2 <= 3
Generals
 v1  v2 
End
