#  Installing IMG on Ubuntu

Installing the GAP IMG package of Laurent Bartholdi on a moderately
recent Ubuntu, starting from the standard desktop installation. The
following assumes that you have administrator rights, in particular,
that you can use the `sudo` command.

Some websites:

-   Ubuntu Linux: <https://ubuntu.com/>
-   Laurent Bartholdi homepage: <https://www.uni-math.gwdg.de/laurent/>
-   GAP (Groups, Algorithms, Programming): <https://www.gap-system.org/>
-   GAP packages:
    -   FGA (Free Group Algorithms):
        <http://www.icm.tu-bs.de/ag_algebra/software/FGA/>
    -   FR (Functionally Recursive Groups):
        <https://gap-packages.github.io/fr/>
    -   GBNP (Gröbner bases of noncommutative polynomials)
        <http://dam03.win.tue.nl/gbnp/chap0.html> -- Download link
		currently broken
	-   IMG (Iterated Monodromy Groups):
        <https://gap-packages.github.io/img/>
    -   LPRES (Nilpotent Quotients of L-Presented Groups)
        <https://gap-packages.github.io/lpres/>
    -   NQ (Nilpotent Quotients of Finitely Presented Groups)
        <https://gap-packages.github.io/nq/>

GAP is included in Ubuntu, but not installed by default. Start by
installing it, together with the recommended packages:

	sudo apt install gap gap-float

Next are some general system and development tools:

	sudo apt install g++ cmake autoconf nodejs

We need a few additional libraries:

	sudo apt install libreadline-dev zlib1g-dev libsuitesparse-dev

You can of course combine these all into one `sudo apt install`
command, or install it via your favorite installation tool.

Now it is time to install local GAP packages. The recommended procedure
given here installs them in your home directory, so other users will not
be able to use these. If you want to install these systemwide, you need
to install these in some global directory for GAP packages. By default,
local GAP packages reside in `$HOME/.gap/pkg/`, and if this directory
does not exist, you need to create it:

	mkdir -p $HOME/.gap/pkg/

Here and below, you can replace `$HOME` by `~`, but there are different
tilde characters, creating problems when copy-pasting.

Now we need to download and unpack the GAP packages FGA, FR, GBNP, IMG,
LPRES, and NQ, which are not included with Ubuntu. (Currently omitting
GBNP because the link seems to be broken.) The last two commands
implicitly assume that you have no other files ending in .tar.gz in the
same directory. To be safe, you can run all these commands in some
temporary directory with nothing important in it.


	wget http://www.icm.tu-bs.de/ag_algebra/software/FGA/FGA-1.4.0.tar.gz
	wget https://github.com/gap-packages/fr/releases/download/v2.4.8/fr-2.4.8.tar.gz
	wget https://github.com/gap-packages/img/releases/download/v0.3.2/IMG-0.3.2.tar.gz
	wget https://github.com/gap-packages/lpres/releases/download/v1.0.1/lpres-1.0.1.tar.gz
	wget https://github.com/gap-packages/nq/releases/download/v2.5.6/nq-2.5.6.tar.gz
	for a in *.tar.gz; do tar xvzf $a -C $HOME/.gap/pkg/ && rm $a; done
	rm *.tar.gz

Obviously, if `wget` fails, you will have to try to get this package in
some other way.

Both the NQ and the IMG package have some code that needs to be
compiled. First IMG:


	cd $HOME/.gap/pkg/IMG-0.3.2/
	aclocal && autoconf && automake
	./configure --with-gaproot=/usr/lib/gap/
	make

Now NQ:

	cd $HOME/.gap/pkg/nq-2.5.6/
	aclocal && autoconf && automake
	./configure --with-gaproot=/usr/lib/gap/
	make

If everything succeeded, you should now be ready to use the IMG
package. (In the following code block, the `$` is the shell command
prompt, and `gap>` is the command prompt for `gap`. Neither of them
should be typed. Make sure not to forget the semicolon at the end of
each `gap` command.) 

	$ gap
	gap> LoadPackage("IMG");
	gap> dendrite := PolynomialSphereMachine(2,[1/6]);
	gap> P1MapBySphereMachine(dendrite);

This last command should output `<z^2+(0.+1.i_z)>`, denoting the
polynomial *f(z) = z<sup>2</sup>+i*

---

*Lukas Geyer <geyer@montana.edu>*
