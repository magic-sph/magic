Contributing to the code
########################

Checking the consistency of the code
====================================


Modifying/contributing to the code
==================================

* Before commiting your modifications **always** make sure that the auto-tests pass correctly.

* Try to follow the same coding style rules as in the rest of the code:

 1. **Never** use TABS but always SPACES instead
 2. Use 3 spaces for indentation
 3. Never use capital letters for variable declaration
 4. Never use :code:`dimension(len)` for declaring array but rather :code:`real(cp) :: data(len)`
 5. Always use the default precisions when introducing new variables :code:`(cp)`


More on that topic `here <http://www.fortran90.org/src/best-practices.html>`_.
