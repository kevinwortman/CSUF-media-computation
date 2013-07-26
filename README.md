CSUF-media-computation
======================

C++ media computation functionality for the CS0 course taught at
California State University, Fullerton.

This module contains object-oriented classes Color, Image, and Clip
which define a visual color, bitmap image, and audio clip,
respectively.  There are functions to read and write Image objects in
the PPM file format, and read and write Clip objects in WAV format.

There is also a suite of functions that operate on a static Image
instance called a "canvas."  These functions allow programmers to
manipulate images using only int and double data types, loops, and
predefined functions.

Likewise, there are a suite of sound_... functions that operate on a
static Clip object.

The PPM format is supported by standard image viewers and GIMP on
Unix, and by IrfanView on Windows.  WAV is perhaps the most widely
supported audio format.

This code is intended for pedagogical use and so I have prioritized
clarity and consistency over efficiency.  In particular, the Color
class stores color components with the double type, and an Image is a
matrix of Color instances stored in STL containers.  This approach is
consistent with what's usually taught in introductory courses, though
is space-inefficient compared to the packed 24-bit RGB layout that is
typically used in production code.

All code is declared inside the csuf namespace.

The module has no external dependencies beyond the standard C++
library.

The unit tests depend on:
  - Boost Unit Test Framework
  - SCons
  - ImageMagick
  - sox
Run the unit tests with the command
  $ scons test

This is open source software distributed under the terms of the
MIT/X11 license.  See LICENSE for copyright and license information.

The elephant.jpg image used by the unit tests is distributed under the
terms of the Creative Commons Attribution 2.0 Generic (CC-BY-2.0)
license.  It was created by flickr user http2007 and available at
http://www.flickr.com/photos/http2007/1149137981/ .
