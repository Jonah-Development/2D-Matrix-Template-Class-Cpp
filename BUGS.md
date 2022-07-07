# Known Bugs
Here is a list of all currently known bugs:
* ...

# Fixed Bugs

* *Bug: When using Print() or cout() on a matrix of type 'char' or 'unsigned char', the value of '0' wouldn't be printed.* **Fixed:** now you can specify a casting type for printing values. If you have a matrix type 'char' or 'unsigned char', use: myMatrix.cout\<int\>();
