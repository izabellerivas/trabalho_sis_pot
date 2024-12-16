   c   l   o   s   e       a   l   l   ;      
   c   l   e   a   r       a   l   l   ;      
   c   l   c   ;      
      
   %   P   A   S   S   O       1   :       C   O   N   V   E   R   T   E   R       T   U   D   O       P   A   R   A       P   U      
      
   S   b   a   s   e   =   1   0   0   ;      
   S   G   2   _   M   V   A   =   5   0   +   3   0   i   ;                                                                   %       G   e   r   a   �   �   o       n   a       b   a   r   r   a       2       e   m       M   V   A      
   S   G   3   _   M   V   A   =   0   ;                                                                                                       %       G   e   r   a   �   �   o       n   a       b   a   r   r   a       3       e   m       M   V   A      
   S   L   2   _   M   V   A   =   3   0   5   .   6   +   1   4   0   .   2   i   ;                               %       C   a   r   g   a       n   a       b   a   r   r   a       2       e   m       M   V   A      
   S   L   3   _   M   V   A   =   1   3   8   .   6   +   4   5   .   2   i   ;                                       %       C   a   r   g   a       n   a       b   a   r   r   a       3       e   m       M   V   A      
      
   S   G   2   =   S   G   2   _   M   V   A   /   S   b   a   s   e   ;                   %       G   e   r   a   �   �   o       n   a       b   a   r   r   a       2       e   m       p   u   ;      
   S   G   3   =   S   G   3   _   M   V   A   /   S   b   a   s   e   ;                   %       G   e   r   a   �   �   o       n   a       b   a   r   r   a       3       e   m       p   u   ;      
   S   L   2   =   S   L   2   _   M   V   A   /   S   b   a   s   e   ;                           %       C   a   r   g   a       n   a       b   a   r   r   a       2       e   m       p   u   ;      
   S   L   3   =   S   L   3   _   M   V   A   /   S   b   a   s   e   ;                           %       C   a   r   g   a       n   a       b   a   r   r   a       3       e   m       p   u   ;      
      
   S   2   i   n   j   =   S   G   2   -   S   L   2   ;      
   S   3   i   n   j   =   S   G   3   -   S   L   3   ;      
      
   %   D   E   F   I   N   I   �   �   O       D   A       M   A   T   R   I   Z       D   E       A   D   M   I   T   �   N   C   I   A      
   y   1   0   =   0   ;      
   y   2   0   =   0   ;      
   y   3   0   =   0   ;      
   y   1   2   =   1   0   -   2   0   i   ;      
   y   1   3   =   1   0   -   3   0   i   ;      
   y   2   3   =   1   6   -   3   2   i   ;      
   y   2   1   =   y   1   2   ;      
   y   3   1   =   y   1   3   ;      
   y   3   2   =   y   2   3   ;      
      
   Y   1   1   =   y   1   0   +   y   1   2   +   y   1   3   ;      
   Y   2   2   =   y   2   0   +   y   2   1   +   y   2   3   ;      
   Y   3   3   =   y   3   0   +   y   3   1   +   y   3   2   ;      
      
   Y   1   2   =   -   y   1   2   ;      
   Y   1   3   =   -   y   1   3   ;      
   Y   2   1   =   Y   1   2   ;      
   Y   2   3   =   -   y   2   3   ;      
   Y   3   1   =   Y   1   3   ;      
   Y   3   2   =   Y   2   3   ;      
      
   Y   =   [   Y   1   1       Y   1   2       Y   1   3      
                           Y   2   1       Y   2   2       Y   2   3      
                           Y   3   1       Y   3   2       Y   3   3   ]   ;      
      
   %   P   A   S   S   O       3   :       C   �   L   C   U   L   O       I   T   E   R   A   T   I   V   O       D   A   S       T   E   N   S   �   E   S       N   A   S       B   A   R   R   A   S      
   V   1   _   i   n   i   =   1   .   0   5   ;      
   V   2   _   i   n   i   =   1   ;      
   V   3   _   i   n   i   =   1   ;      
      
      
      
   f   o   r       k   =   1   :   5      
   V   1   (   k   )   =   V   1   _   i   n   i   (   k   )   ;      
   V   2   (   k   )   =   V   2   _   i   n   i   (   k   )   ;      
   V   3   (   k   )   =   V   3   _   i   n   i   (   k   )   ;      
      
   V   1   _   i   n   i   (   k   +   1   )   =   V   1   (   k   )   ;      
   V   2   _   i   n   i   (   k   +   1   )   =   (   c   o   n   j   (   S   2   i   n   j   /   V   2   (   k   )   )   -   Y   2   1   *   V   1   (   k   )   -   Y   2   3   *   V   3   (   k   )   )   .   /   Y   2   2   ;      
   V   3   _   i   n   i   (   k   +   1   )   =   (   c   o   n   j   (   S   3   i   n   j   /   V   3   (   k   )   )   -   Y   3   1   *   V   1   (   k   )   -   Y   3   2   *   V   2   _   i   n   i   (   k   +   1   )   )   .   /   Y   3   3   ;      
   e   n   d      
   %   C   �   L   C   U   L   O       D   A       P   O   T   E   N   C   I   A       I   N   J   E   T   A   D   A       N   A       B   A   R   R   A       S   L   A   C   K      
   V   =   [   V   1   '       V   2   '       V   3   '   ]      
   m   o   d   u   l   V   =   a   b   s   (   [   V   1   '       V   2   '       V   3   '   ]   )      
   a   n   g   l   e   V   =   a   n   g   l   e   (   [   V   1   '       V   2   '       V   3   '   ]   )   *   1   8   0   /   p   i      
   V   e   n   d   =   [   V   1   (   e   n   d   )       V   2   (   e   n   d   )       V   3   (   e   n   d   )   ]      
      
      
   V   =   [   V   1   (   e   n   d   )       V   2   (   e   n   d   )       V   3   (   e   n   d   )   ]   ;      
      
      
   i   =   1   ;      
   P   =   0   ;      
   Q   =   0   ;      
   f   o   r       j   =   1   :   3      
   P   =   P   +   a   b   s   (   V   (   i   ,   i   )   )   *   a   b   s   (   V   (   i   ,   j   )   )   *   a   b   s   (   Y   (   i   ,   j   )   )   *   c   o   s   (   a   n   g   l   e   (   Y   (   i   ,   j   )   )   -   a   n   g   l   e   (   V   (   i   ,   i   )   )   +   a   n   g   l   e   (   V   (   i   ,   j   )   )   )   ;      
   Q   =   Q   -   a   b   s   (   V   (   i   ,   i   )   )   *   a   b   s   (   V   (   i   ,   j   )   )   *   a   b   s   (   Y   (   i   ,   j   )   )   *   s   i   n   (   a   n   g   l   e   (   Y   (   i   ,   j   )   )   -   a   n   g   l   e   (   V   (   i   ,   i   )   )   +   a   n   g   l   e   (   V   (   i   ,   j   )   )   )   ;      
   e   n   d      
   S   1   =   [   P       Q   ]   ;      
      
   %   C   A   L   C   U   L   O       D   A   S       P   E   R   D   A   S       N   A   S       L   I   N   H   A   S      
   P   p   =   z   e   r   o   s   (   3   ,   3   )      
      
   f   o   r       i   =   1   :   3      
                   f   o   r       j   =   1   :   3      
                   P   p   (   i   ,   j   )   =   -   a   b   s   (   V   (   i   )   )   ^   2   *   a   b   s   (   Y   (   i   ,   j   )   )   *   c   o   s   (   a   n   g   l   e   (   Y   (   i   ,   j   )   )   )   +   a   b   s   (   V   (   i   )   )   *   a   b   s   (   V   (   j   )   )   *   a   b   s   (   Y   (   i   ,   j   )   )   *   c   o   s   (   a   n   g   l   e   (   Y   (   i   ,   j   )   )   -   a   n   g   l   e   (   V   (   i   )   )   +   a   n   g   l   e   (   V   (   j   )   )   )   ;      
                   Q   p   (   i   ,   j   )   =   a   b   s   (   V   (   i   )   )   ^   2   *   a   b   s   (   Y   (   i   ,   j   )   )   *   s   i   n   (   a   n   g   l   e   (   Y   (   i   ,   j   )   )   )   -   a   b   s   (   V   (   i   )   )   *   a   b   s   (   V   (   j   )   )   *   a   b   s   (   Y   (   i   ,   j   )   )   *   s   i   n   (   a   n   g   l   e   (   Y   (   i   ,   j   )   )   -   a   n   g   l   e   (   V   (   i   )   )   +   a   n   g   l   e   (   V   (   j   )   )   )   ;      
                   i   f       j   =   =   i      
                                   P   p   (   i   ,   j   )   =   0   ;      
                   e   n   d      
                   e   n   d      
   e   n   d      
      
      
   P   p      
   Q   p      
      
   P   p   l   i   n   h   a   s   =   [   P   p   (   1   ,   2   )   +   P   p   (   2   ,   1   )      
                                           P   p   (   1   ,   3   )   +   P   p   (   3   ,   1   )      
                                           P   p   (   2   ,   3   )   +   P   p   (   3   ,   2   )   ]      
      
   Q   p   l   i   n   h   a   s   =   [   Q   p   (   1   ,   2   )   +   Q   p   (   2   ,   1   )      
                                           Q   p   (   1   ,   3   )   +   Q   p   (   3   ,   1   )      
                                           Q   p   (   2   ,   3   )   +   Q   p   (   3   ,   2   )   ]      
