.. include:: <isonum.txt>

============
Dependencies
============

STaGE is released under GNU LGPL V3.0 and may be released under a later version of
the LGPL. The external programs used by STaGE are released under their own licences.
Some of them cannot be distributed with STaGE, but must be downloaded on their own.

`ACPYPE <https://code.google.com/p/acpype/>`_ (required)
    Used for generating GAFF and OPLS-AA output.

    * `Sousa da Silva, A. W., Vranken, W. ACPYPE - AnteChamber PYthon Parser interfacE.
      BMC Res. Notes 2012, 5: 367.
      <http://www.biomedcentral.com/1756-0500/5/367>`_

`AnteChamber <http://ambermd.org/antechamber/antechamber.html>`_ (required)
    Used for generating GAFF and OPLS-AA output.

    * `Wang, J., Wang, W., Kollman P. A., Case, D. A. Automatic atom type and
      bond type perception in molecular mechanical calculations.
      J. Mol. Graph. and Model. 2006, 25: 247-260.
      <http://www.sciencedirect.com.ezp.sub.su.se/science/article/pii/S1093326305001737>`_
    * `Wang, J., Wolf, R. M., Caldwell, J. W., Kollman, P. A., Case, D. A.
      Development and testing of a general AMBER force field.
      J. Comput. Chem. 2004, 25: 1157-1174.
      <http://onlinelibrary.wiley.com.ezp.sub.su.se/doi/10.1002/jcc.20035/abstract>`_

`GAMESS/US <http://www.msg.ameslab.gov/gamess/>`_
    Used for B3LYP/PCM charge calculations.

    * `Schmidt, M. W., Baldridge, K. K., Boatz, J. A., Elbert, S. T., Gordon,
      M. S., Jensen, J. H., Koseki, S., Matsunaga, N., Nguyen, K. A., Su, S.,
      Windus, T. L., Dupuis, M., Montgomery, J. A. General atomic and
      molecular electronic structure system. J. Comput. Chem. 1993,
      14: 1347-1363. <http://onlinelibrary.wiley.com.ezp.sub.su.se/doi/10.1002/jcc.540141112/abstract>`_

gmstoresp.sh
    Included in STaGE distribution. Used for B3LYP/PCM charge calculations.

    |copy| 2004, Sarnoff Corporation, Princeton, NJ, USA

`GROMACS <http://www.gromacs.org/>`_
    Used for solvating and neutralizing (if needed) the system.

    * `Pronk, S., Páll, S., Schulz, R., Larsson, P., Bjelkmar, P., Apostolov, R.,
      Shirts, M. R., Smith, J. C., Kasson, P. M., van der Spoel, D., Hess, B., Lindahl, E.
      GROMACS 4.5: a high-throughput and highly parallel open source molecular
      simulation toolkit. Bioinformatics 2013, 29: 845-854.
      <http://bioinformatics.oxfordjournals.org.ezp.sub.su.se/content/29/7/845>`_
    * `Hess, B., Kutzner, C., van der Spoel, D., Lindahl, E. GROMACS 4: Algorithms
      for Highly Efficient, Load-Balanced, and Scalable Molecular Simulation.
      J. Chem. Theory Comput. 2008, 4: 435-447.
      <http://pubs.acs.org.ezp.sub.su.se/doi/abs/10.1021/ct700301q>`_

`MATCH <http://brooks.chem.lsa.umich.edu/index.php?page=match&subdir=articles/resources/software>`_
    Used for generating CGenFF output.

    * `Yesselman, J. D., Price, D. J., Knight, J. L. Brooks, C. L. 3rd. MATCH: an
      atom-typing toolset for molecular mechanics force fields.
      J. Comput. Chem. 2012, 33: 189-202.
      <http://onlinelibrary.wiley.com.ezp.sub.su.se/doi/10.1002/jcc.21963/abstract>`_

`Open Babel <http://openbabel.org/>`_ (required)
    Used for molecule format conversions and for assigning OPLS-AA atom types.

    * `O' Boyle, M. O., Banck, M., James, C. A., Morley, C., Vandermeersch, T.,
      Hutchison, G. R. Open Babel: An open checmical toolbox.
      J. Cheminform. 2011, 3: 33.
      <http://www.jcheminf.com/content/3/1/33>`_

---------------------
Environment Variables
---------------------

The following applications must be available in the list of directories of
executable programs (the PATH environment variable):

    * acpype.py
    * antechamber
    * amsol7.1
    * tleap
    * babel
    * MATCH.pl
    * pdb2gmx
    * editconf
    * genbox
    * grompp
    * genion
    * make_ndx
    * genrestr
    * rungms
    * gmstoresp.sh

    The following environment variables must be set

    * AMBERHOME
    * MATCH
    * PerlChemistry
