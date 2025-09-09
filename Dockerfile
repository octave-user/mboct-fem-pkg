# syntax=docker/dockerfile:1

## Copyright (C) 2024(-2024) Reinhard <octave-user@a1.net>

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

#####################################################################
## ABOUT THIRD PARTY SOFTWARE
#####################################################################

## Any docker image which is built from this Dockerfile may contain
## additional third party software in binary form. All third party
## software is distributed under their respective license and terms
## of use.

## Some of those copyright statements are included below.

## Despite of extensive research by the author,
## there is no guarantee that this section is really comprehensive.
## So, just in case that the reader may find any third-party software
## inside such an image, but the respective license statement is not yet
## included here, please inform the author about the license text
## to be included!

#####################################################################
## GNU-Octave
#####################################################################

## Octave is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <https://www.gnu.org/licenses/>.

#####################################################################
## TFEL
#####################################################################

## Description

## `TFEL` is a collaborative development of CEA and EDF.

## `MFront` is a code generator which translates a set of closely related
## domain specific languages into plain C++ on top of the `TFEL`
## library. Those languages covers three kind of material knowledge:

## - material properties (for instance the
##   Young modulus, the thermal conductivity, etc.)
## - mechanical behaviours. Numerical performances of generated
##   mechanical behaviours was given a particular attention. MFront
##   offers a variety of interfaces to finite element solvers `Cast3M`,
##   `Code-Aster`, `EUROPLEXUS`, `Abaqus-Standard`, `Abaqus-Explicit`,
##   `Zebulon`, etc.. or various FFT solvers such as
##   `AMITEX_FFTP`. Various benchmarks shows that `MFront`
##   implementations are competitive with native implementations
##   available in the `Cast3M`, `Code-Aster` and `Cyrano3` solvers.
## - simple point-wise models, such as material swelling
##   used in fuel performance codes.

## `MFront` comes with an handy easy-to-use tool called `MTest` that can
## test the local behaviour of a material, by imposing independent
## constraints on each component of the strain or the stress. This tool
## has been much faster (from ten to several hundred times depending on
## the test case) than using a full-fledged finite element solver.

## # Licences

## `TFEL` version prior to 0.1 were released under both the LGPL and the
## CECILL-B licences. A copy of those licences are included in the
## distributions of TFEL.

## `TFEL` versions 1.x were developped by CEA within the PLEIADES
## project. Since svn revision 584, TFEL was part of the `PLEIADES`
## projet.

## Starting from versions 2.x, TFEL has been publicly released under
## either the GPL or the CECILL-A licence. A copy of thoses licences are
## delivered with the sources of TFEL. CEA or EDF may also distribute
## this project under specific licensing conditions.

## Copyright (C) 2006-2013 CEA/DEN. All rights reserved.
## Copyright (C) 2014-2015 CEA/DEN, EDF R&D. All rights reserved.

#####################################################################
## MFrontGenericInterfaceSupport
#####################################################################

## This project aims at providing tools (functions, classes, bindings,
## etc...) to handle behaviours written using `MFront` generic interface.
## For information about `MFront`, see
## <http://thelfer.github.io/tfel/web/index.html>.

## Those tools are meant to be used by (`FEM`, `FFT`, etc.) solver
## developers. This tools are *not* linked to the `TFEL` libraries.
## Permissive licences have been chosen to allow integration in open-source
## and proprietary codes.

## This project is described in this paper:
## [![DOI](https://joss.theoj.org/papers/10.21105/joss.02003/status.svg)](https://doi.org/10.21105/joss.02003)

## The official website can be found here:
## <https://thelfer.github.io/mgis/web/index.html>.

## ## The `MFrontGenericInterface` `C++` library

## The project is build around the `MFrontGenericInterface` library. This
## library provides two main functions:

## - the `mgis::behaviour::load` functions loads `MFront` behaviours from
##   external shared libraries and retrieve all relevant meta data
##   function. Those relevant information are stored in the
##   `mgis::behaviour::Behaviour` class.
## - the `mgis::behaviour::integrate` integrates the behaviour over one
##   time step. The data associated with an integration point are handled
##   by the `mgis::behaviour::BehaviourData` class which contains the state
##   of the integration point at the beginning and at the end of the time
##   step.

## The library also supports handling a group of integration points though
## the `mgis::behaviour::MaterialStateManager` class.

## ## Bindings

## ### Existing

## The following bindings are available:

## - `c` bindings.
## - `python` bindings .
## - `fortran03` bindings.
## - `fenics` bindings (under current work in the `master` branch). Those
##   bindings are strongly inspired by the `fenics-solid-mechanics`
##   project. Those bindings are currently quite limited as mostly serve
##   as a proof of concept. Note that `MGIS` can also be used in `FEniCS`
##   through the `python` interface. This is discussed here:
##   <https://thelfer.github.io/mgis/web/FEniCSBindings.html>.
## - `julia` bindings (experimental)

## ### Future bindings (contributors are welcomed)

## The following bindings are under consideration:

## - `octave` binding

## # Versions, branches

## - Version `2.2` is meant to be build against `TFEL` 4.2
## - Version `2.1` is meant to be build against `TFEL` 4.1
## - Version `2.0` is meant to be build against `TFEL` 4.0
## - Version `1.2.2` is meant to be build against `TFEL` 3.4.3
## - Version `1.2.1` is meant to be build against `TFEL` 3.4.1
## - Version `1.2` is meant to be build against `TFEL` 3.4.0
## - Version `1.1` is meant to be build against `TFEL` 3.3.0
## - Version `1.0` is meant to be build against `TFEL` 3.2.0
## - Version `1.0.1` is meant to be build against `TFEL` 3.2.1

## The following branches are available:

## - The `master` branch follows the evolution of the `master` branch of
##   the `TFEL` project
## - The `rliv-2.2` follows the evolution of the 4.2.x series of the `TFEL`
##   project.
## - The `rliv-2.1` follows the evolution of the 4.1.x series of the `TFEL`
##   project.
## - The `rliv-2.0` follows the evolution of the 4.0.x series of the `TFEL`
##   project.
## - The `rliv-1.2` follows the evolution of the 3.4.x series of the `TFEL`
##   project.
## - The `rliv-1.1` follows the evolution of the 3.3.x series of the `TFEL`
##   project.
## - The `rliv-1.0` follows the evolution of the 3.2.x series of the `TFEL`
##   project. Note that this branch is **not** compatible with
##   `TFEL-3.2.0`.

## # Acknowledgement

## This project uses code extracted from the following projects:

## - https://github.com/bitwizeshift/string_view-standalone by Matthew
##   Rodusek
## - https://github.com/mpark/variant: by Michael Park
## - https://github.com/progschj/ThreadPool by Jakob Progsch and Vaclav
##   Zeman
## - https://github.com/martinmoene/span-lite by Martin Moene
## - https://bitbucket.org/fenics-apps/fenics-solid-mechanics/ by
##   Kristian B. Olgaard and Garth N. Wells.

##  Use, modification and distribution are subject
##  to one of the following licences:
##  - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
##    file LGPL-3.0.txt)
##  - CECILL-C,  Version 1.0 (See accompanying files
##    CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).

#####################################################################
## MFrontGallery
#####################################################################

## Contributing to the `MFrontGallery` project

## Contributions to `MFrontGallery` are greatly appreciated.

## Please take a moment to review this document in order to make the
## contribution process easy and effective for everyone involved.

## Following these guidelines helps to communicate that you respect the time of
## the developers managing and developing this open source project. In return,
## they should reciprocate that respect in addressing your issue or assessing
## patches and features.

## ## Using the issue tracker

## The [issue
## tracker](https://github.com/thelfer/MFrontGallery/issues)
## is the preferred channel for [bug reports](#bugs), [features
## requests](#features) and [submitting pull requests](#pull-requests), but
## please respect the following restrictions:

## * Please **do not** use the issue tracker for personal support requests
##   (contact directly the authors
##   [tfel-contact@cea.fr](mailto:tfel-contact@cea.fr)

## ## Bug reports

## A bug is a _demonstrable problem_ that is caused by the code in the repository.
## Good bug reports are extremely helpful - thank you!

## Guidelines for bug reports:

## 1. **Use the GitHub issue search**: check if the issue has already been
##    reported.

## 2. **Check if the issue has been fixed**: try to reproduce it using the
##    latest `master` or development branch in the repository.

## 3. **Isolate the problem**: ideally create a [reduced test
##    case].

## A good bug report shouldn't leave others needing to chase you up for more
## information. Please try to be as detailed as possible in your report. What is
## your environment? What steps will reproduce the issue? What compiler(s) and OS
## experience the problem? What would you expect to be the outcome? All these
## details will help people to fix any potential bugs.

## Example:

## > Short and descriptive example bug report title
## >
## > A summary of the issue, versions of `TFEL/MFront` used and the
## > OS/compiler environment in which it occurs. If suitable, include the
## > steps required to reproduce the bug.
## >
## > 1. This is the first step
## > 2. This is the second step
## > 3. Further steps, etc.
## >
## > Any other information you want to share that is relevant to the issue being
## > reported. This might include the lines of code that you have identified as
## > causing the bug, and potential solutions (and your opinions on their
## > merits).

## ## Feature requests

## Feature requests are welcome. But take a moment to find out whether your idea
## fits with the scope and aims of the project. It's up to *you* to make a strong
## case to convince the project's developers of the merits of this feature. Please
## provide as much detail and context as possible.

## ## Pull requests

## Good pull requests - patches, improvements, new features - are a fantastic
## help. They should remain focused in scope and avoid containing unrelated
## commits.

## **Please ask first** before embarking on any significant pull request (e.g.
## implementing features, refactoring code, porting to a different language),
## otherwise you risk spending a lot of time working on something that the
## project's developers might not want to merge into the project.

## Please adhere to the coding conventions used throughout a project (indentation,
## accurate comments, etc.) and any other requirements (such as test coverage).

## Adhering to the following this process is the best way to get your work
## included in the project:

## 1. [Fork](http://help.github.com/fork-a-repo/) the project, clone your fork,
##    and configure the remotes:

##    ```bash
##    # Clone your fork of the repo into the current directory
##    git clone https://github.com/<your-username>/MFrontGallery.git
##    # Navigate to the newly cloned directory
##    cd MFrontGallery
##    # Assign the original repo to a remote called "upstream"
##    git remote add upstream https://github.com/thelfer/MFrontGallery.git
##    ```

## 2. If you cloned a while ago, get the latest changes from upstream:

##    ```bash
##    git checkout master
##    git pull upstream master
##    ```

## 3. Create a new topic branch (off the main project development branch) to
##    contain your feature, change, or fix:

##    ```bash
##    git checkout -b <topic-branch-name>
##    ```

## 4. Commit your changes in logical chunks. Please adhere to these [git commit
##    message guidelines](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
##    or your code is unlikely be merged into the main project. Use Git's
##    [interactive rebase](https://help.github.com/articles/interactive-rebase)
##    feature to tidy up your commits before making them public.

## 5. Locally merge (or rebase) the upstream development branch into your topic branch:

##    ```bash
##    git pull [--rebase] upstream master
##    ```

## 6. Push your topic branch up to your fork:

##    ```bash
##    git push origin <topic-branch-name>
##    ```

## 7. [Open a Pull Request](https://help.github.com/articles/using-pull-requests/)
##     with a clear title and description.

## **IMPORTANT**: By submitting a patch, you agree to allow the project owners to
## license your work under the the terms of the *LGPL License*.

## <!--
## This file is merely a copy of the `CONTRIBUTING.md` file of the
## html5boilerplate project

## Copyright (c) HTML5 Boilerplate

## Permission is hereby granted, free of charge, to any person obtaining a copy of
## this software and associated documentation files (the "Software"), to deal in
## the Software without restriction, including without limitation the rights to
## use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
## of the Software, and to permit persons to whom the Software is furnished to do
## so, subject to the following conditions:

## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
## -->

#####################################################################
## NLopt
#####################################################################

## NLopt combines several free/open-source nonlinear optimization
## libraries by various authors.  See the COPYING, COPYRIGHT, and README
## files in the subdirectories for the original copyright and licensing
## information of these packages.

## The compiled NLopt library, i.e. the combined work of all of the
## included optimization routines, is licensed under the conjunction of
## all of these licensing terms.  Currently, the most restrictive terms
## are for the code in the "luksan" directory, which is licensed under
## the GNU Lesser General Public License (GNU LGPL), version 2.1 or
## later (see luksan/COPYRIGHT).

## That means that the compiled NLopt library is governed by the terms of
## the LGPL.

## ---------------------------------------------------------------------------

## Other portions of NLopt, including any modifications to the abovementioned
## packages, are licensed under the standard "MIT License:"

## Copyright (c) 2007-2011 Massachusetts Institute of Technology

## Permission is hereby granted, free of charge, to any person obtaining
## a copy of this software and associated documentation files (the
## "Software"), to deal in the Software without restriction, including
## without limitation the rights to use, copy, modify, merge, publish,
## distribute, sublicense, and/or sell copies of the Software, and to
## permit persons to whom the Software is furnished to do so, subject to
## the following conditions:

## The above copyright notice and this permission notice shall be
## included in all copies or substantial portions of the Software.

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
## LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
## WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#####################################################################
## Trilinos
#####################################################################

## ************************************************************************
##
##            Trilinos: An Object-Oriented Solver Framework
##                 Copyright (2001) Sandia Corporation
##
##
## Copyright (2001) Sandia Corporation. Under the terms of Contract
## DE-AC04-94AL85000, there is a non-exclusive license for use of this
## work by or on behalf of the U.S. Government.  Export of this program
## may require a license from the United States Government.
##
## 1. Redistributions of source code must retain the above copyright
## notice, this list of conditions and the following disclaimer.
##
## 2. Redistributions in binary form must reproduce the above copyright
## notice, this list of conditions and the following disclaimer in the
## documentation and/or other materials provided with the distribution.
##
## 3. Neither the name of the Corporation nor the names of the
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
## EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
## CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
## EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
## PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
## PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
## NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## NOTICE:  The United States Government is granted for itself and others
## acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
## license in this data to reproduce, prepare derivative works, and
## perform publicly and display publicly.  Beginning five (5) years from
## July 25, 2001, the United States Government is granted for itself and
## others acting on its behalf a paid-up, nonexclusive, irrevocable
## worldwide license in this data to reproduce, prepare derivative works,
## distribute copies to the public, perform publicly and display
## publicly, and to permit others to do so.
##
## NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
## OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
## ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
## RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
## INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
## THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
##

#####################################################################
## googletest
#####################################################################

## Copyright 2008, Google Inc.
## All rights reserved.

## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:

##     * Redistributions of source code must retain the above copyright
## notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above
## copyright notice, this list of conditions and the following disclaimer
## in the documentation and/or other materials provided with the
## distribution.
##     * Neither the name of Google Inc. nor the names of its
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.

## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#####################################################################
## MBDyn
#####################################################################

## MBDyn (C) is a multibody analysis code.
## http://www.mbdyn.org
##
## Copyright (C) 1996-2023
##
## Pierangelo Masarati     <pierangelo.masarati@polimi.it>
## Paolo Mantegazza        <paolo.mantegazza@polimi.it>
##
## Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
## via La Masa, 34 - 20156 Milano, Italy
## http://www.aero.polimi.it
##
## Changing this copyright notice is forbidden.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation (version 2 of the License).
##
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##
## ----------------------------------------------------------------

#####################################################################
## NetCDF
#####################################################################

## Copyright 2018 Unidata

## Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

## 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

## 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

## 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## The NetCDF Copyright.

## \page copyright Copyright

## Copyright 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
## 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 University
## Corporation for Atmospheric Research/Unidata.

## Portions of this software were developed by the Unidata Program at the
## University Corporation for Atmospheric Research.

## Access and use of this software shall impose the following obligations
## and understandings on the user. The user is granted the right, without
## any fee or cost, to use, copy, modify, alter, enhance and distribute
## this software, and any derivative works thereof, and its supporting
## documentation for any purpose whatsoever, provided that this entire
## notice appears in all copies of the software, derivative works and
## supporting documentation.  Further, UCAR requests that the user credit
## UCAR/Unidata in any publications that result from the use of this
## software or in any product that includes this software, although this
## is not an obligation. The names UCAR and/or Unidata, however, may not
## be used in any advertising or publicity to endorse or promote any
## products or commercial entity unless specific written permission is
## obtained from UCAR/Unidata. The user also understands that
## UCAR/Unidata is not obligated to provide the user with any support,
## consulting, training or assistance of any kind with regard to the use,
## operation and performance of this software nor to provide the user
## with any updates, revisions, new versions or "bug fixes."

## THIS SOFTWARE IS PROVIDED BY UCAR/UNIDATA "AS IS" AND ANY EXPRESS OR
## IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL UCAR/UNIDATA BE LIABLE FOR ANY SPECIAL,
## INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
## FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
## NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
## WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.


#####################################################################
## Intel MKL
#####################################################################

## Intel Simplified Software License (Version October 2022)
## Use and Redistribution. You may use and redistribute the software, which is
## provided in binary form only, (the "Software"), without modification, provided the
## following conditions are met:
## * Redistributions must reproduce the above copyright notice and these terms of use
## in the Software and in the documentation and/or other materials provided with
## the distribution.
## * Neither the name of Intel nor the names of its suppliers may be used to endorse
## or promote products derived from this Software without specific prior written
## permission.
## * No reverse engineering, decompilation, or disassembly of the Software is
## permitted, nor any modification or alteration of the Software or its operation
## at any time, including during execution.
## No other licenses. Except as provided in the preceding section, Intel grants no
## licenses or other rights by implication, estoppel or otherwise to, patent,
## copyright, trademark, trade name, service mark or other intellectual property
## licenses or rights of Intel.
## Third party software. "Third Party Software" means the files (if any) listed in
## the "third-party-software.txt" or other similarly-named text file that may be
## included with the Software. Third Party Software, even if included with the
## distribution of the Software, may be governed by separate license terms, including
## without limitation, third party license terms, open source software notices and
## terms, and/or other Intel software license terms. These separate license terms
## solely govern Your use of the Third Party Software.
## DISCLAIMER. THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
## WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT ARE
## DISCLAIMED. THIS SOFTWARE IS NOT INTENDED FOR USE IN SYSTEMS OR APPLICATIONS WHERE
## FAILURE OF THE SOFTWARE MAY CAUSE PERSONAL INJURY OR DEATH AND YOU AGREE THAT YOU
## ARE FULLY RESPONSIBLE FOR ANY CLAIMS, COSTS, DAMAGES, EXPENSES, AND ATTORNEYS'
## FEES ARISING OUT OF ANY SUCH USE, EVEN IF ANY CLAIM ALLEGES THAT INTEL WAS
## NEGLIGENT REGARDING THE DESIGN OR MANUFACTURE OF THE SOFTWARE.
## LIMITATION OF LIABILITY. IN NO EVENT WILL INTEL BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
## NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
## PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
## WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.
## No support. Intel may make changes to the Software, at any time without notice,
## and is not obligated to support, update or provide training for the Software.
## Termination. Your right to use the Software is terminated in the event of your
## breach of this license.
## Feedback. Should you provide Intel with comments, modifications, corrections,
## enhancements or other input ("Feedback") related to the Software, Intel will be
## free to use, disclose, reproduce, license or otherwise distribute or exploit the
## Feedback in its sole discretion without any obligations or restrictions of any
## kind, including without limitation, intellectual property rights or licensing
## obligations.
## Compliance with laws. You agree to comply with all relevant laws and regulations
## governing your use, transfer, import or export (or prohibition thereof) of the
## Software.
## Governing law. All disputes will be governed by the laws of the United States of
## America and the State of Delaware without reference to conflict of law principles
## and subject to the exclusive jurisdiction of the state or federal courts sitting
## in the State of Delaware, and each party agrees that it submits to the personal
## jurisdiction and venue of those courts and waives any objections. THE UNITED
## NATIONS CONVENTION ON CONTRACTS FOR THE INTERNATIONAL SALE OF GOODS (1980) IS
## SPECIFICALLY EXCLUDED AND WILL NOT APPLY TO THE SOFTWARE.

FROM ubuntu:latest

LABEL org.opencontainers.image.title="mboct-fem-pkg"
LABEL org.opencontainers.image.vendor="Reinhard"
LABEL com.docker.desktop.extension.icon="https://yt3.googleusercontent.com/HNK0rqhf4dCLnlP99g8c1odImPBZLM_QxLRb7yGBGDlDCbiiReQmoRZK4fuSU0XLC4yGApVS5Q=s160-c-k-c0x00ffffff-no-rj"
LABEL org.opencontainers.image.description="This package belongs to a suite of packages which can be used for pre- and postprocessing of flexible bodies in MBDyn (www.mbdyn.org) with GNU-Octave. It contains a general purpose structural Finite Element toolkit for linear statics and dynamics, which can be used to generate flexible body data for MBDyn's modal element and also hydrodynamic lubricated journal and slider plain bearings."
LABEL com.docker.extension.screenshots=[{"alt"="Fourbar","url"="https://i.ytimg.com/an_webp/d4i5AYPxsG4/mqdefault_6s.webp?du=3000&sqp=CI6PubUG&rs=AOn4CLDf2JH21dC1U_B8UytYHSNBP_LV9g"},{"alt"="Twisted bar","url"="https://i.ytimg.com/an_webp/I8HENx5mszA/mqdefault_6s.webp?du=3000&sqp=CICNubUG&rs=AOn4CLCr0TUUL4UZoynMQxoO7JMz27qNRg"},{"alt"="Cook's membrane","url"="https://i.ytimg.com/vi/rxQP8V4U0dE/hqdefault.jpg?sqp=-oaymwE2CPYBEIoBSFXyq4qpAygIARUAAIhCGAFwAcABBvABAfgBsASAAuADigIMCAAQARgTIB4ofzAP&rs=AOn4CLC7YAoIVIleUQWKwJ00c3tNb-1CBg"}]
LABEL org.opencontainers.image.source=https://github.com/octave-user
LABEL com.docker.extension.publisher-url="https://github.com/octave-user"
LABEL com.docker.extension.additional-urls=[{"title":"MBDyn","url":"https://www.mbdyn.org"},{"title":"MBDyn-GitLab","url":"https://public.gitlab.polimi.it/DAER/mbdyn"},{"title":"Fourbar","url":"https://www.youtube.com/watch?v=d4i5AYPxsG4"},{"title":"Twisted bar","url":"https://www.youtube.com/watch?v=I8HENx5mszA"},{"title":"Rolling ring","url":"https://www.youtube.com/watch?v=rxQP8V4U0dE"},{"title":"Cook's membrane","url":"https://www.youtube.com/watch?v=EAgejp4jQ00"},{"title":"videos","url":"https://www.youtube.com/@nonlinearmultibodydynamics6802"}]
LABEL org.opencontainers.image.licenses=GPL3

ARG SRC_DIR=/usr/local/src/
ARG LICENSE_DIR=/usr/local/share/license/
ARG BUILD_DIR=/tmp/build/
ARG TESTS_DIR=/tmp/tests/
ARG MBD_NUM_TASKS=4
ARG RUN_TESTS='octave;trilinos;mbdyn;mboct'
ARG RUN_CONFIGURE='none'
ARG CXX=g++
ARG CC=gcc
ARG FC=gfortran
ARG F77=gfortran
ARG AWKPATH="${AWKPATH}:${BUILD_DIR}"

WORKDIR ${SRC_DIR}
WORKDIR ${BUILD_DIR}

COPY Dockerfile ${SRC_DIR}
COPY Dockerfile ${BUILD_DIR}
COPY octave-source.awk ${BUILD_DIR}
COPY octave-source.sh ${BUILD_DIR}

RUN sed 's/Types: deb/Types: deb deb-src/g' -i /etc/apt/sources.list.d/ubuntu.sources

RUN rm -f /etc/apt/apt.conf.d/docker-clean; echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' > /etc/apt/apt.conf.d/keep-cache

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get -yq update && apt-get -yq build-dep octave && \
    apt-get -yq install git libopenmpi-dev libmumps-seq-dev \
    libnlopt-dev libhdf5-dev libginac-dev wget libmetis-dev \
    libnetcdf-c++4-dev parallel cmake clang++-19

ARG GMSH_URL="http://www.gmsh.info/bin/Linux/"
ARG GMSH_VERSION="stable"
ENV GMSH_TAR="gmsh-${GMSH_VERSION}-Linux64.tgz"

WORKDIR ${SRC_DIR}/gmsh
WORKDIR ${BUILD_DIR}/gmsh

RUN --mount=type=cache,target=${BUILD_DIR}/gmsh,sharing=locked <<EOT bash
    if ! test -f "${GMSH_TAR}"; then
      if ! wget "${GMSH_URL}${GMSH_TAR}"; then
        exit 1
      fi
    fi

    if ! test -d gmsh-*.*.*-Linux64/bin; then
      if ! tar -zxvf "${GMSH_TAR}"; then
        exit 1
      fi
    fi

    if ! install gmsh-*.*.*-Linux64/bin/gmsh /usr/local/bin; then
      exit 1
    fi

    if ! gmsh --version; then
      exit 1
    fi
    cp "${GMSH_TAR}" "${SRC_DIR}/gmsh"
EOT

WORKDIR ${SRC_DIR}/tfel
WORKDIR ${BUILD_DIR}/tfel

ARG TFEL_REPO=https://github.com/thelfer/tfel.git
ENV LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

RUN --mount=type=cache,target=${BUILD_DIR}/tfel,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/tfel/.git; then
      if ! git clone ${TFEL_REPO} ${BUILD_DIR}/tfel; then
        exit 1
      fi
    fi

    if ! git pull --force; then
      exit 1
    fi

    if ! test -d build_dir; then
      mkdir build_dir
    fi

    cd build_dir

    if ! test -f Makefile; then
      if ! cmake -S .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-19 -DCMAKE_C_COMPILER=clang-19; then
        exit 1
      fi
    fi

    if ! make -j${MBD_NUM_TASKS}; then
      exit 1
    fi

    case "${RUN_TESTS}" in
    *tfel*|all)
      if ! make -j${MBD_NUM_TASKS} check; then
        exit 1
      fi
      ;;
    none)
      ;;
    esac

    if ! make install; then
      exit 1
    fi

    make package_source

    find . -name '*-Source.tar.gz' -exec cp '{}' ${SRC_DIR}/tfel ';'

    if test -f install_manifest.txt; then
      cat install_manifest.txt | awk '/\/lib.*\.so$/' | xargs chmod +x
    fi
EOT

ARG MGIS_REPO=https://github.com/thelfer/MFrontGenericInterfaceSupport.git

WORKDIR ${SRC_DIR}/mgis
WORKDIR ${BUILD_DIR}/mgis

RUN --mount=type=cache,target=${BUILD_DIR}/mgis,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/mgis/.git; then
      if ! git clone ${MGIS_REPO} ${BUILD_DIR}/mgis; then
        exit 1
      fi
    fi

    if ! git pull --force; then
      exit 1
    fi

    if ! test -d build_dir; then
      mkdir build_dir
    fi

    cd build_dir

    if ! test -f Makefile; then
      if ! cmake -DTFELTests_DIR=/usr/local/share/tfel/cmake -DTFELException_DIR=/usr/local/share/tfel/cmake -DTFELUtilities_DIR=/usr/local/share/tfel/cmake -DTFELMaterial_DIR=/usr/local/share/tfel/cmake -DMTestFileGenerator_DIR=/usr/local/share/tfel/cmake -S ..; then
        exit 1
      fi
    fi

    if ! make -j${MBD_NUM_TASKS} all; then
      exit 1
    fi

    case "${RUN_TESTS}" in
      *mgis*|all)
      if ! make -j${MBD_NUM_TASKS} check; then
        exit 1
      fi
      ;;
    none)
      ;;
    esac

    if ! make install; then
      exit 1
    fi

    make package_source

    find . -name '*-Source.tar.gz' -exec cp '{}' ${SRC_DIR}/mgis ';'

    if test -f install_manifest.txt; then
      cat install_manifest.txt | awk '/\/lib.*\.so$/' | xargs chmod +x
    fi
EOT

WORKDIR ${SRC_DIR}/gallery
WORKDIR ${BUILD_DIR}/gallery

ARG GALLERY_REPO="https://github.com/thelfer/MFrontGallery.git"
ARG GALLERY_BRANCH="master"
RUN --mount=type=cache,target=${BUILD_DIR}/gallery,sharing=locked <<EOT bash
    if ! test -d "${BUILD_DIR}/gallery/.git"; then
      if ! git clone -b "${GALLERY_BRANCH}" "${GALLERY_REPO}" "${BUILD_DIR}/gallery"; then
        exit 1
      fi
    fi

    if ! git pull --force; then
      exit 1
    fi

    if ! test -d build_dir; then
      mkdir build_dir
    fi

    cd build_dir

    if ! test -f Makefile; then
      if ! cmake -S .. -DCMAKE_BUILD_TYPE=Release -Denable-generic=ON -Denable-generic-behaviours=ON; then
        echo "Warning: cmake failed"
      fi
    fi

    if ! make -j${MBD_NUM_TASKS} all; then
      exit 1
    fi

    case "${RUN_TESTS}" in
      *gallery*|all)
      if ! make -j${MBD_NUM_TASKS} check; then
        exit 1
      fi
      ;;
    none)
      ;;
    esac

    if ! make install; then
      exit 1
    fi

    make package_source

    find . -name '*-Source.tar.gz' -exec cp '{}' ${SRC_DIR}/gallery ';'

    if test -f install_manifest.txt; then
      cat install_manifest.txt | awk '/\/lib.*\.so$/' | xargs chmod +x
    fi
EOT

WORKDIR ${SRC_DIR}/octave
WORKDIR ${BUILD_DIR}/octave

RUN --mount=type=cache,target=${BUILD_DIR}/octave,sharing=locked <<EOT bash
    "${BUILD_DIR}/octave-source.sh"
EOT

WORKDIR ${LICENSE_DIR}/mkl
RUN wget --output-document=mkl-license.pdf https://cdrdv2.intel.com/v1/dl/getContent/749362

WORKDIR ${BUILD_DIR}

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get -yq install libmkl-full-dev

# WORKDIR ${SRC_DIR}/siconos
# WORKDIR ${BUILD_DIR}/siconos

# ENV SICONOS_REPOSITORY="https://github.com/siconos/siconos.git"
# ENV SICONOS_BRANCH="master"

# RUN --mount=type=cache,target=${BUILD_DIR}/siconos,sharing=locked <<EOT bash
#     if ! test -d ${BUILD_DIR}/siconos/.git; then
#       git clone -b "${SICONOS_BRANCH}" "${SICONOS_REPOSITORY}" "${BUILD_DIR}/siconos"
#     fi

#     git pull --force

#     if ! test -d build_dir; then
#       mkdir build_dir
#     fi

#     cd build_dir

#     if ! test -f Makefile; then
#       export CPPFLAGS="`python3-config --cflags` ${CPPFLAGS}"
#       export LDFLAGS="`python3-config --ldflags` ${LDFLAGS}"
#       export LIBS="`python3-config --libs` ${LIBS}"
#       cmake -S .. -DSICONOS_SYSTEM_WIDE_INSTALL=ON -DWITH_TESTING=OFF  -DWITH_PYTHON_WRAPPER=OFF -DWITH_CMAKE_BUILD_TYPE=Release ## -DUSER_OPTIONS_FILE=../config_samples/siconos_lite.cmake
#     fi

#     make -j${MBD_NUM_TASKS} all
#     make install
# EOT

WORKDIR ${SRC_DIR}/Trilinos
WORKDIR ${BUILD_DIR}/Trilinos

ARG TRILINOS_REPO="https://github.com/trilinos/Trilinos.git"
ARG TRILINOS_BRANCH="master"
ARG TRILINOS_CONFIG="-DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DTrilinos_ENABLE_NOX=ON -DTrilinos_ENABLE_Epetra=ON -DTrilinos_ENABLE_EpetraExt=ON -DTrilinos_ENABLE_Amesos=ON -DTrilinos_ENABLE_AztecOO=ON -DEpetra_ENABLE_MPI=OFF -DNOX_ENABLE_Epetra=ON -DNOX_ENABLE_EpetraExt=ON -DNOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA=ON -DNOX_ENABLE_AztecOO=ON -DNOX_ENABLE_Ifpack=ON -DTrilinos_ENABLE_TESTS=OFF"
ARG TRILINOS_PREFIX="/usr/local/"

RUN --mount=type=cache,target=${BUILD_DIR}/Trilinos,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/Trilinos/.git; then
      git clone -b ${TRILINOS_BRANCH} ${TRILINOS_REPO} ${BUILD_DIR}/Trilinos
    fi

    cd ${BUILD_DIR}/Trilinos

    if ! test -d build_dir; then
      mkdir build_dir
    fi

    cd build_dir

    case "${RUN_CONFIGURE}" in
    *trilinos*|all)
      rm -f Makefile
      ;;
    none)
      ;;
    esac

    if ! test -f Makefile; then
      cmake .. -DCMAKE_INSTALL_PREFIX="${TRILINOS_PREFIX}" ${TRILINOS_CONFIG}
    fi

    make -j${MBD_NUM_TASKS}

    case "${RUN_TESTS}" in
      *trilinos*|all)
      if ! make -j${MBD_NUM_TASKS} check; then
        exit 1
      fi
      ;;
    none)
      ;;
    esac

    ## FIXME: make package_source requires too much disc space
    # make package_source

    # find . -name '*-Source.tar.gz' -exec cp '{}' ${SRC_DIR}/Trilinos ';'

    make install
EOT

# WORKDIR ${SRC_DIR}/pastix
# WORKDIR ${BUILD_DIR}/pastix

# ARG PASTIX_REPO="https://gitlab.inria.fr/solverstack/pastix.git"
# ARG PASTIX_BRANCH="master"
# ARG PASTIX_CONFIG="-DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release"
# ARG PASTIX_PREFIX="/usr/local/"

# RUN --mount=type=cache,target=${BUILD_DIR}/pastix,sharing=locked <<EOT bash
#     if ! test -d ${BUILD_DIR}/pastix/.git; then
#       git clone -b ${PASTIX_BRANCH} ${PASTIX_REPO} ${BUILD_DIR}/pastix
#     fi

#     cd ${BUILD_DIR}/pastix

#     if ! test -d build_dir; then
#       mkdir build_dir
#     fi

#     cd build_dir

#     case "${RUN_CONFIGURE}" in
#     *pastix*|all)
#       rm -f Makefile
#       ;;
#     none)
#       ;;
#     esac

#     if ! test -f Makefile; then
#       cmake .. -DCMAKE_INSTALL_PREFIX="${PASTIX_PREFIX}" ${PASTIX_CONFIG}
#     fi

#     make -j${MBD_NUM_TASKS}

#     case "${RUN_TESTS}" in
#       *pastix*|all)
#       if ! make -j${MBD_NUM_TASKS} check; then
#         exit 1
#       fi
#       ;;
#     none)
#       ;;
#     esac

#     make package_source

#     find . -name '*-Source.tar.gz' -exec cp '{}' ${SRC_DIR}/pastix ';'

#     make install
# EOT

WORKDIR ${SRC_DIR}/gtest
WORKDIR ${BUILD_DIR}/gtest

ARG GTEST_REPO="https://github.com/google/googletest.git"
ARG GTEST_BRANCH="main"

RUN --mount=type=cache,target=${BUILD_DIR}/gtest,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/gtest/.git; then
      git clone -b ${GTEST_BRANCH} ${GTEST_REPO} ${BUILD_DIR}/gtest
    fi

    cd ${BUILD_DIR}/gtest

    if ! test -d build_dir; then
      mkdir build_dir
    fi

    cd build_dir

    case "${RUN_CONFIGURE}" in
    *gtest*|all)
      rm -f Makefile
      ;;
    none)
      ;;
    esac

    if ! test -f Makefile; then
      cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
    fi

    make -j${MBD_NUM_TASKS}

    make package_source

    find . -name '*-Source.tar.gz' -exec cp '{}' ${SRC_DIR}/gtest ';'

    make install
EOT

ARG MBD_FLAGS="-Ofast -Wall -march=native -mtune=native"
ARG MBD_CPPFLAGS="-I/usr/local/include -I/usr/local/include/kokkos -I/usr/include/suitesparse -I/usr/include/mkl -I/usr/local/include/MGIS -I/usr/local/include/MFront"
ARG MBD_CXXFLAGS=""
ARG MBD_ARGS_WITH="--with-mfront --with-static-modules --with-arpack --with-umfpack --with-klu --with-arpack --with-lapack --without-metis --without-mpi --with-trilinos --with-pardiso --with-suitesparseqr --with-qrupdate --with-gtest"
ARG MBD_ARGS_ENABLE="--enable-octave --enable-multithread --disable-Werror --enable-install_test_progs"
ARG MBD_REPO="https://public.gitlab.polimi.it/DAER/mbdyn.git"

WORKDIR ${SRC_DIR}/mbdyn
WORKDIR ${BUILD_DIR}/mbdyn

RUN --mount=type=cache,target=${BUILD_DIR}/mbdyn,sharing=locked <<EOT bash
    if ! test -d ${BUILD_DIR}/mbdyn/.git; then
      if ! git clone -b develop ${MBD_REPO} ${BUILD_DIR}/mbdyn; then
        exit 1
      fi
    fi

    git pull --force

    if ! test -x ./configure; then
      ./bootstrap.sh
    fi

    case "${RUN_CONFIGURE}" in
    *mbdyn*|all)
      rm -f Makefile
      ;;
    none)
      ;;
    esac

    if ! test -f Makefile; then
      ./configure CPPFLAGS="${MBD_CPPFLAGS}" LDFLAGS="${MBD_FLAGS}" CXXFLAGS="${MBD_FLAGS} ${MBD_CXXFLAGS}" CFLAGS="${MBD_FLAGS}" FCFLAGS="${MBD_FLAGS}" F77FLAGS="${MBD_FLAGS}" ${MBD_ARGS_WITH} ${MBD_ARGS_ENABLE}
    fi

    if ! make -j${MBD_NUM_TASKS} all; then
      exit 1
    fi

    case "${RUN_TESTS}" in
    *mbdyn*|all)
      if ! make test; then
        exit 1
      fi
      ;;
    none)
      ;;
    esac

    if ! make install; then
      exit 1
    fi

    if ! /usr/local/mbdyn/bin/mbdyn --version; then
      exit 1
    fi

    if make dist; then
      cp mbdyn-*.tar.gz ${SRC_DIR}/mbdyn
    fi
EOT

WORKDIR ${SRC_DIR}/octave-pkg
WORKDIR ${BUILD_DIR}/octave-pkg

ARG MBOCT_OCTAVE_PKG_REPO="https://github.com/octave-user/mboct-octave-pkg.git"
ARG MBOCT_OCTAVE_PKG_BRANCH="master"
ARG MBOCT_NUMERICAL_PKG_REPO="https://github.com/octave-user/mboct-numerical-pkg.git"
ARG MBOCT_NUMERICAL_PKG_BRANCH="master"
ARG MBOCT_MBDYN_PKG_REPO="https://github.com/octave-user/mboct-mbdyn-pkg.git"
ARG MBOCT_MBDYN_PKG_BRANCH="master"
ARG MBOCT_FEM_PKG_REPO="https://github.com/octave-user/mboct-fem-pkg.git"
ARG MBOCT_FEM_PKG_BRANCH="master"

RUN --mount=type=cache,target=${BUILD_DIR}/octave-pkg,sharing=locked <<EOT bash
    octave -q --eval 'pkg install -verbose -global -forge nurbs;pkg install -verbose -global -forge netcdf'

    if ! test -d mboct-octave-pkg; then
      git clone -b "${MBOCT_OCTAVE_PKG_BRANCH}" "${MBOCT_OCTAVE_PKG_REPO}" 'mboct-octave-pkg'
    fi

    pushd mboct-octave-pkg && git pull --force && popd

    make CXXFLAGS="${MBD_FLAGS}" -j${MBD_NUM_TASKS} -C 'mboct-octave-pkg' dist install_global

    if ! test -d mboct-numerical-pkg; then
      git clone -b "${MBOCT_NUMERICAL_PKG_BRANCH}" "${MBOCT_NUMERICAL_PKG_REPO}" mboct-numerical-pkg
    fi

    pushd mboct-numerical-pkg && git pull --force && popd

    make CXXFLAGS="${MBD_FLAGS}" -j${MBD_NUM_TASKS} -C 'mboct-numerical-pkg' dist install_global

    if ! test -d mboct-mbdyn-pkg; then
      git clone -b "${MBOCT_MBDYN_PKG_BRANCH}" "${MBOCT_MBDYN_PKG_REPO}" mboct-mbdyn-pkg
    fi

    pushd mboct-mbdyn-pkg && git pull --force && popd

    make CXXFLAGS="${MBD_FLAGS}" -j${MBD_NUM_TASKS} -C 'mboct-mbdyn-pkg' dist install_global

    if ! test -d mboct-fem-pkg; then
      git clone -b "${MBOCT_FEM_PKG_BRANCH}" "${MBOCT_FEM_PKG_REPO}" mboct-fem-pkg
    fi

    pushd mboct-fem-pkg && git pull --force && popd

    make CXXFLAGS="${MBD_FLAGS}" -j${MBD_NUM_TASKS} -C 'mboct-fem-pkg' dist install_global

    pushd mboct-fem-pkg/src

    if ! test -x configure; then
      ./bootstrap
    fi

    case "${RUN_CONFIGURE}" in
    *mboct*|all)
      rm -f Makefile
      ;;
    none)
      ;;
    esac

    if ! test -f Makefile; then
      ./configure CXXFLAGS="${MBD_FLAGS}"
    fi

    make -j${MBD_NUM_TASKS} all
    make install
    make distclean

    popd

    find . -name 'mboct-*-pkg-*.tar.gz' -exec cp '{}' ${SRC_DIR}/octave-pkg ';'
EOT

WORKDIR ${BUILD_DIR}/octave-pkg

RUN --mount=type=cache,target=${BUILD_DIR}/octave-pkg,sharing=locked <<EOT bash
    case "${RUN_TESTS}" in
    *mboct*|all)
      make MBD_NUM_TASKS=2 -C mboct-octave-pkg check_installed_parallel
      make MBD_NUM_TASKS=2 -C mboct-numerical-pkg check_installed_parallel
      make MBD_NUM_TASKS=2 -C mboct-mbdyn-pkg check_installed_parallel
      make MBD_NUM_TASKS=2 -C mboct-fem-pkg check_installed_parallel
      ;;
    none)
      ;;
    esac
EOT

WORKDIR /home/ubuntu
## RUN rm -rf ${BUILD_DIR} ${TESTS_DIR} ## Clean up temporary files
USER ubuntu
ENV LANG=C

## Run on Ubuntu with graphics enabled
## xhost + local:docker; docker run --entrypoint octave --user root -it --network=host --env="HOME=$HOME" --workdir="$HOME" --volume="/run/user:/run/user:rw" --volume="$HOME:$HOME:rw" --rm -v /tmp/.X11-unix:/tmp/.X11-unix:rw -e DISPLAY=${DISPLAY} --env="XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR}" --env=LIBGL_ALWAYS_INDIRECT=0 --env=LIBGL_ALWAYS_SOFTWARE=1 -h $HOSTNAME -v $HOME/.Xauthority:/home/ubuntu/.Xauthority octaveuser/mboct-fem-pkg:latest --persist --norc --eval 'pkg("local_list","/home/ubuntu/octave_packages");'; xhost -;

ENTRYPOINT [ "/usr/local/bin/octave", "--persist", "--norc", "--eval", "pkg('local_list','/home/ubuntu/octave_packages');pkg('load','mboct-fem-pkg');" ]
