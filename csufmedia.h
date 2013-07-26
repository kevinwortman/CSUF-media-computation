
//////////////////////////////////////////////////////////////////////
// csufmedia.h
//
// C++ media computation functionality for the CS0 course taught at
// California State University, Fullerton.
//
// See LICENSE.txt for copyright and license information.
//
// Coding Conventions ////
//
// A color intensity represents a red, green, or blue intensity and
// should be between 0.0 and 1.0 inclusive.
//
// Dimensions follow the standard graphics convention:
// - the top-left corner is (0,0)
// - the bottom-right corner is (width-1, height-1)
//
// (0, 0) +----+
//        |    |
//        +----+ (width-1, height-1)
//
// Thus valid x coordinates are between 0 and width-1, inclusive; and
// valid y coordinates are between 0 and height-1, inclusive.
//
// An audio sample represents an amplitude value between -1.0 and +1.0
// inclusive.
//
// Function arguments are checked with assert() calls.
//
// I/O functions return true on success and false on failure.
//
// The code loosely conforms to the Google C++ Style Guide.
//
// References ////
//
// http://en.wikipedia.org/wiki/Netpbm_format
// https://ccrma.stanford.edu/courses/422/projects/WaveFormat/
// http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml
// 
//////////////////////////////////////////////////////////////////////

#ifndef _CSUFMEDIA_H
#define _CSUFMEDIA_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace csuf {

  //////////////////////////////////////////////////////////////////////
  // Constants
  //////////////////////////////////////////////////////////////////////

  // sample rate used by the sound_ functions
  const int SOUND_SAMPLE_RATE = 44100;

  //////////////////////////////////////////////////////////////////////
  // Inline utility functions
  //////////////////////////////////////////////////////////////////////

  // Return true if and only if x is between least and most,
  // inclusive.
  inline bool is_in_range(int least, int x, int most) {
    return (least <= x) && (x <= most);
  }
  inline bool is_in_range(double least, double x, double most) {
    return (least <= x) && (x <= most);
  }

  // Return true if and only if x is a valid color intensity,
  // i.e. between 0.0 and 1.0, inclusive.
  inline bool is_valid_intensity(double x) {
    return is_in_range(0.0, x, 1.0);
  }

  // Return true if and only if x is a valid audio sample,
  // i.e. between -1.0 and +1.0 inclusive.
  inline bool is_valid_sample(double x) {
    return is_in_range(-1.0, x, 1.0);
  }

  // Return true if an only if x is a valid delta value, between 0.0
  // and 1.0 inclusive.
  inline bool is_valid_delta(double x) {
    return is_valid_intensity(x);
  }

  //////////////////////////////////////////////////////////////////////
  // Canvas functions
  //////////////////////////////////////////////////////////////////////

  // Set the canvas to be a blank image with the given width and
  // height, which must both be positive.  Every pixel in the canvas
  // is initialized to black.
  void make_canvas(int width, int height);

  // Load an image in the PPM format into the canvas.  Any prior
  // contents of the canvas are overwritten.  The function returns
  // true on success and false on error.  If an error is encountered,
  // the canvas is left empty.
  bool load_canvas(const std::string& ppm_file_path);

  // Save the canvas to a file in the PPM format.  The canvas must not
  // be empty.
  bool save_canvas(const std::string& ppm_file_path);

  // Return true when the canvas is empty, false otherwise.
  bool is_canvas_empty();

  // Get the canvas dimensions.  If the canvas is empty, 0 is
  // returned.
  int canvas_width();
  int canvas_height();
  
  // Get a color intensity value at the given (x, y) coordinate.
  double get_red(int x, int y);
  double get_green(int x, int y);
  double get_blue(int x, int y);

  // Set a color intensity value at the given (x, y) coordinate.
  double set_red(int x, int y, double red);
  double set_green(int x, int y, double green);
  double set_blue(int x, int y, double blue);

  // Set all three color intensities at the given (x, y) coordinates.
  void set_pixel(int x, int y, double red, double green, double blue);

  // Draw a filled rectangle.
  void rectangle(int x, int y,
		 int width, int height,
		 double red, double green, double blue);

  // Fill a vertical semicircle.  More precisely, 180 degrees of the
  // area between two concentric circles (an "annulus").  (x, y) is
  // the top-left corner of the rectangular region covered by the
  // semicircle.  If rightward is true then the circle center is
  // (x, y+height/2); otherwise it is (x+width, y+height/2).
  void vertical_semicircle(int x, int y,
			   int inner_radius, int outer_radius,
			   bool rightward,
			   double red, double green, double blue);

  //////////////////////////////////////////////////////////////////////
  // Sound functions
  //////////////////////////////////////////////////////////////////////

  // Create a sound with the given number of samples.  All samples are
  // initialized to 0.0 (silence).
  void make_sound(int sample_count);

  // Load the current sound from the given WAV file.  If the file has
  // more than one channel, they are downmixed into one channel.
  // Returns true on success.
  bool load_sound(const std::string& wav_file_path);

  // Save the current sound to the given WAV file.  The sound must not
  // be empty.  Returns true on success.
  bool save_sound(const std::string& wav_file_path);

  // Returns true if the current sound is empty.
  bool is_sound_empty();

  // Returns the number of samples in the current sound, which must
  // not be empty.
  int sound_sample_count();

  // Get one sample from the current sound, which must not be empty.
  double get_sample(int i);

  // Set one sample in the current sound, which must not be empty.
  void set_sample(int i, double sample);

  bool repeat_sound(double samples[],
		    int sample_count,
		    int repeat_count,
		    const std::string& wav_file_path);

  //////////////////////////////////////////////////////////////////////
  // Object-oriented image interface
  //////////////////////////////////////////////////////////////////////

  ////
  // A Color represents a visible color as a tuple of red, green, blue
  // intensity values.
  ////
  class Color {
  public:
    Color(double red, double green, double blue)
      : _red(red),
	_green(green),
	_blue(blue) {
      assert(is_valid_intensity(red));
      assert(is_valid_intensity(green));
      assert(is_valid_intensity(blue));
    }

    ////
    // Operators
    ////

    bool operator== (const Color& o) const {
      return ((_red == o._red) &&
	      (_green == o._green) &&
	      (_blue == o._blue));
    }

    ////
    // Accessors
    ////

    double get_red()   const { return _red;   }
    double get_green() const { return _green; }
    double get_blue()  const { return _blue;  }

    ////
    // Mutators
    ////

    void set_red(double red) {
      assert(is_valid_intensity(red));
      _red = red;
    }

    void set_green(double green) {
      assert(is_valid_intensity(green));
      _green = green;
    }

    void set_blue(double blue) {
      assert(is_valid_intensity(blue));
      _blue = blue;
    }

  private:
    double _red, _green, _blue;
  };

  ////
  // Color constants
  ////

  const Color BLACK(0, 0, 0),
    WHITE(1, 1, 1),
    RED(1, 0, 0),
    GREEN(0, 1, 0),
    BLUE(0, 0, 1);

  ////
  // An Image is an image represented by a rectangular grid of pixels,
  // each stored as a Color object.  An Image may be empty or
  // nonempty; when empty, an Image contains no pixels and few of its
  // methods work.
  ////
  class Image {
  public:
    // Create an empty Image.
    Image() {
      assert(is_empty());
    }

    // Create an Image with the given dimensions, and all pixels
    // initialized to the given color.
    Image(int width, int height, const Color& color) {
      resize(width, height, color);
      assert(!is_empty());
    }

    ////
    // Operators
    ////

    bool operator== (const Image& o) const {
      return _rows == o._rows;
    }

    ////
    // Accessors
    ////

    // Returns true if the Image is empty.
    bool is_empty() const {
      return _rows.empty();
    }

    // Get the dimensions of the Image.  The image must not be empty.
    int get_height() const {
      assert(!is_empty());
      return static_cast<int>(_rows.size());
    }
    int get_width() const {
      assert(!is_empty());
      return static_cast<int>(_rows.front().size());
    }

    // Returns true if the given coordinates are valid.  If the image
    // is empty, all coordinates are invalid.
    bool is_valid_coordinates(int x, int y) const {
      return (!is_empty()
	      && is_in_range(0, x, get_width() - 1)
	      && is_in_range(0, y, get_height() - 1));
    }

    ////
    // Mutators for the entire image
    ////

    // Make the Image empty.
    void clear() {
      _rows.clear();
      assert(is_empty());
    }

    // Change to the given dimensions.  The current contents of the
    // image is replaced with pixels of the given color.
    void resize(int width, int height, const Color& color) {
      assert(width > 0);
      assert(height > 0);
      _rows.assign(height, std::vector<Color>(width, color));
      assert(!is_empty());
    }

    // Replace all pixels with the given color.
    void fill(const Color& color) {
      // Note: it's OK for the image to be empty, in that case the
      // following loop never iterates.
      for (std::vector<std::vector<Color > >::iterator i = _rows.begin();
	   i != _rows.end();
	   ++i)
	i->assign(get_width(), color);
    }
	
    ////
    // Mutators for individual pixels
    ////

    // Get the color of the pixel at the given coordinates.
    const Color& get_pixel(int x, int y) const {
      assert(is_valid_coordinates(x, y));
      return _rows[y][x];
    }

    // Get the color of the pixel at the given coordinates, as a
    // mutable reference which may be modified.
    Color& get_mutable_pixel(int x, int y) {
      assert(is_valid_coordinates(x, y));
      return _rows[y][x];
    }
      
    // Set the color of the pixel at the given coordinates.
    void set_pixel(int x, int y, const Color& color) {
      assert(is_valid_coordinates(x, y));
      _rows[y][x] = color;
    }

  private:
    // 2D matrix of Color objects
    std::vector<std::vector<Color> > _rows;
  };

  //////////////////////////////////////////////////////////////////////
  // Object-oriented audio interface
  //////////////////////////////////////////////////////////////////////

  ////
  // A Clip is a sequence of audio samples.  It may contain multiple
  // channels; 1 channel is monoraul and 2 channels is stereo.  A Clip
  // knows its sample rate, but doesn't know any other metadata.  A
  // Clip object can be empty, in which case few of its methods work.
  ////
  class Clip {
  public:
    // Create an empty clip.
    Clip() {
      assert(is_empty());
    }

    ////
    // Accessors
    ////

    bool is_empty() const {
      return _samples.empty();
    }

    int channel_count() const {
      return static_cast<int>(_samples.size());
    }

    int is_channel_index(int x) const {
      return is_in_range(0, x, channel_count() - 1);
    }

    int frame_count() const {
      if (is_empty())
	return 0;
      else
	return static_cast<int>(_samples.front().size());
    }

    int is_frame_index(int x) const {
      return is_in_range(0, x, frame_count() - 1);
    }

    double get_sample(int channel, int frame) const {
      assert(!is_empty());
      assert(is_channel_index(channel));
      assert(is_frame_index(frame));
      return _samples[channel][frame];
    }

    int get_sample_rate() const {
      assert(!is_empty());
      return _sample_rate;
    }

    ////
    // Mutators for the entire clip
    ////

    void assign(int channel_count, int frame_count, int sample_rate, double default_sample) {
      assert(channel_count > 0);
      assert(frame_count > 0);
      assert(sample_rate > 0);
      assert(is_valid_sample(default_sample));
      
      _samples.assign(channel_count, std::vector<double>(frame_count, default_sample));
      _sample_rate = sample_rate;
      
      assert(!is_empty());
    }

    void clear() {
      _samples.clear();
      assert(is_empty());
    }

    void downmix();
    	
    ////
    // Mutators for individual samples
    ////

    void set_sample(int channel, int frame, double sample) {
      assert(!is_empty());
      assert(is_channel_index(channel));
      assert(is_frame_index(frame));
      assert(is_valid_sample(sample));
      _samples[channel][frame] = sample;
    }

  private:
    std::vector<std::vector<double> > _samples;
    int _sample_rate;
  };

  //////////////////////////////////////////////////////////////////////
  // Rasterization
  //////////////////////////////////////////////////////////////////////

  // Fill a rectangle.
  void rasterize_rectangle(Image& target,
			   int x, int y,
			   int width, int height,
			   const Color& color);

  // Fill a vertical semicircle.  More precisely, 180 degrees of the
  // area between two concentric circles (an "annulus").  (x, y) is
  // the top-left corner of the rectangular region covered by the
  // semicircle.  If rightward is true then the circle center is
  // (x, y+height/2); otherwise it is (x+width, y+height/2).
  void rasterize_vertical_semicircle(Image& target,
				     int x, int y,
				     int inner_radius, int outer_radius,
				     bool rightward,
				     const Color& color);

  //////////////////////////////////////////////////////////////////////
  // Approximate comparison
  //////////////////////////////////////////////////////////////////////

  // Return true if all the intensity values of the two images are
  // within delta of each other.  This is a way of comparing two
  // images that may've been modified by some process involving small
  // multiplicative errors.  For example, converting between file
  // seems to introduce small rounding errors.
  bool approximate_equal(const Image& a, const Image& b, double delta);

  bool approximate_equal(const Clip& a, const Clip& b, double delta);

  //////////////////////////////////////////////////////////////////////
  // PPM File Format
  //////////////////////////////////////////////////////////////////////

  // Write the given Image into the output stream using the PPM
  // format.  If binary is true, write a binary PPM, otherwise write
  // an ASCII PPM.  If verbose is true, write status reports to cerr.
  // Returns true on success and false on error.
  bool write_ppm(const Image& image,
		 std::ostream& out,
		 bool binary,
		 bool verbose);

  // Read an Image in PPM format out of the input stream.  The image
  // may be in either ASCII or binary (raw) PPM format.  If verbose is
  // true, write status reports to cerr.  Returns true on success and
  // false on error.  In the error case the image object is made
  // empty.
  bool read_ppm(Image& image, 
		std::istream& in,
		bool verbose);

  // Write a PPM file at the given path.  Inputs and outputs are the
  // same as write_ppm.
  bool write_ppm_file(const Image& image,
		      const std::string& path,
		      bool binary,
		      bool verbose);

  // Read a PPM file from the given path.  Inputs and outputs are the
  // same as write_ppm.
  bool read_ppm_file(Image& image,
		     const std::string& path,
		     bool verbose);

  //////////////////////////////////////////////////////////////////////
  // WAV File Format
  //////////////////////////////////////////////////////////////////////

  // Write the given clip to the given output stream in the WAV
  // format.  bytes_per_sample represents the byte depth of the
  // samples, and must be 1 or 2, for 8- or 16-bit samples,
  // respectively.  Returns true on success.
  bool write_wav(const Clip& clip,
		 std::ostream& out,
		 int bytes_per_sample,
		 bool verbose);

  // Read a WAV file from the given input stream.  The file must be
  // encoded using only linear PCM quantization (the most common
  // kind).  Returns true on success.
  bool read_wav(Clip& clip,
		std::istream& in,
		bool verbose);

  bool write_wav_file(const Clip& clip,
		      const std::string& path,
		      int bytes_per_sample,
		      bool verbose);

  bool read_wav_file(Clip& clip,
		     const std::string& path,
		     bool verbose);
}

#endif
