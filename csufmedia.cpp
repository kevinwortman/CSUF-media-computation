
//////////////////////////////////////////////////////////////////////
// csufmedia.cpp
//
// See csufmedia.h for substantive comments.
//
// See LICENSE.txt for copyright and license information.
// 
//////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>

#include "csufmedia.h"

using namespace std;

namespace csuf {

  //////////////////////////////////////////////////////////////////////
  // Constants
  //////////////////////////////////////////////////////////////////////

  const string PPM_ASCII_HEADER = "P3",
               PPM_BINARY_HEADER = "P6";

  const string WAV_CHUNK_ID = "RIFF",
    WAV_FORMAT = "WAVE",
    WAV_SUB_CHUNK_1_ID = "fmt ",
    WAV_SUB_CHUNK_2_ID = "data";
  const int WAV_AUDIO_FORMAT = 1;

  const int PPM_COLUMN_LIMIT = 70;

  const int SOUND_BYTES_PER_SAMPLE = 2;

  //////////////////////////////////////////////////////////////////////
  // Static variables
  //////////////////////////////////////////////////////////////////////

  static Image _canvas;
  static Clip _sound;

  //////////////////////////////////////////////////////////////////////
  // Inline utility functions
  //////////////////////////////////////////////////////////////////////

  inline bool approx_equal(double a, double b, double delta) {
    assert(is_valid_intensity(a) || is_valid_sample(a));
    assert(is_valid_intensity(b) || is_valid_sample(b));
    assert(is_valid_delta(delta));
    return fabs(a - b) <= delta;
  }

  inline bool approx_equal(const Color& a, const Color& b, double delta) {
    return (approx_equal(a.get_red(), b.get_red(), delta) &&
	    approx_equal(a.get_green(), b.get_green(), delta) &&
	    approx_equal(a.get_blue(), b.get_blue(), delta));
  }

  inline bool write_uint(ostream& f, unsigned x, int byte_count) {
    assert(byte_count >= 1);
    assert(byte_count <= 4);

    for (int i = 0; i < byte_count; ++i) {
      f.put((x >> (i*8)) & 0xFF);
      if (f.fail())
	return false;
    }

    return true;
  }

  inline bool read_uint(istream& f, unsigned& x, int byte_count) {
    assert(byte_count >= 1);
    assert(byte_count <= 4);

    x = 0;
    for (int i = 0; i < byte_count; ++i) {
      unsigned byte = f.get();
      if (f.fail())
	return false;
      x |= byte << (i*8);
    }

    return true;
  }

  // Write a string to the given stream, without using the stream
  // insertion operator.
  inline bool write_string(ostream& f, const string& s) {
    f.write(s.c_str(), s.size());
    return f.good();
  }

  // Write a color intensity to a PPM file in either binary or ASCII
  // format.
  inline void write_intensity(ostream& out,
			      int& column,
			      double intensity,
			      bool binary) {
    assert(column >= 0);
    assert(column <= PPM_COLUMN_LIMIT);
    assert(is_valid_intensity(intensity));

    int out_of_255 = static_cast<int>(intensity * 255.0);
    assert(out_of_255 >= 0);
    assert(out_of_255 <= 255);

    if (binary) {
      out.put(out_of_255);
    } else {
      stringstream ss(stringstream::out);
      ss << setw(4) << out_of_255;
      string s = ss.str();

      if ((column + static_cast<int>(s.size())) > PPM_COLUMN_LIMIT) {
	out.put('\n');
	column = 0;
      }

      write_string(out, s);
      column += s.size();
    }
  }

  // Read the next integer from the given input stream of ASCII
  // characters.  Leading whitespace is skipped.  Comments starting
  // with # are skipped to the end of the line.  Returns true on
  // success and false on error.
  inline bool read_int(istream& f, int& x) {
    // first skip anything prior to digits
    while (true) {
      char ch = f.peek();
      if (!f.good())
	return false; // I/O error
      else if (isspace(ch))
	f.get(); // skip whitespace
      else if (isdigit(ch))
	break; // found digits
      else if (ch == '#') {
	// skip to end of line
	while (f.get() != '\n')
	  if (!f.good())
	    return false;
      } else
	return false; // parse error
    }

    // read digits in least-significant-first order
    x = 0;
    while (isdigit(f.peek())) {
      if (!f.good()) // I/O error
	return false;
      else
	x = x*10 + (f.get() - '0');
    }
    
    return true;
  }

  inline bool read_fixed_string(istream& f, string& s, int length) {
    assert(length > 0);
    s.resize(length);
    for (int i = 0; i < length; ++i)
      s[i] = f.get();
    return f.good();
  }

  // Convert an int in the range [0, saturated] to a double in the range
  // [0, 1.0].
  inline double byte_to_intensity(int x, int saturated) {
    assert(x >= 0);
    assert(x <= saturated);
    assert(saturated > 0);
    return static_cast<double>(x) / static_cast<double>(saturated);
  }

  inline string str(int x) {
    ostringstream s;
    s << x;
    return s.str();
  }

  inline void log(bool verbose,
		  const std::string& prefix,
		  const std::string& message) {
    if (verbose)
      cerr << prefix << ": " << message << endl;
  }

  inline void log_ppm(bool verbose,
		      const std::string& message) {
    log(verbose, "ppm", message);
  }

  inline void log_wav(bool verbose,
		      const std::string& message) {
    log(verbose, "wav", message);
  }

  //////////////////////////////////////////////////////////////////////
  // Canvas functions
  //////////////////////////////////////////////////////////////////////

  void make_canvas(int width, int height) {
    _canvas = Image(width, height, BLACK);
  }

  bool load_canvas(const string& ppm_file_path) {
    return read_ppm_file(_canvas, ppm_file_path, true);
  }

  bool save_canvas(const string& ppm_file_path) {
    assert(!is_canvas_empty());
    return write_ppm_file(_canvas, ppm_file_path, true, true);
  }

  bool is_canvas_empty() {
    return _canvas.is_empty();
  }

  int canvas_width() {
    assert(!is_canvas_empty());
    return _canvas.get_width();
  }

  int canvas_height() {
    assert(!is_canvas_empty());
    return _canvas.get_height();
  }

  double get_red(int x, int y) {
    assert(!is_canvas_empty());
    return _canvas.get_pixel(x, y).get_red();
  }

  double get_green(int x, int y) {
    assert(!is_canvas_empty());
    return _canvas.get_pixel(x, y).get_green();
  }

  double get_blue(int x, int y) {
    assert(!is_canvas_empty());
    return _canvas.get_pixel(x, y).get_blue();
  }

  double set_red(int x, int y, double red) {
    assert(!is_canvas_empty());
    _canvas.get_mutable_pixel(x, y).set_red(red);
    return red;
  }

  double set_green(int x, int y, double green) {
    assert(!is_canvas_empty());
    _canvas.get_mutable_pixel(x, y).set_green(green);
    return green;
  }

  double set_blue(int x, int y, double blue) {
    assert(!is_canvas_empty());
    _canvas.get_mutable_pixel(x, y).set_blue(blue);
    return blue;
  }

  void set_pixel(int x, int y, double red, double green, double blue) {
    assert(!is_canvas_empty());
    Color& p = _canvas.get_mutable_pixel(x, y);
    p.set_red(red);
    p.set_green(green);
    p.set_blue(blue);
  }

  void rectangle(int x, int y,
		 int width, int height,
		 double red, double green, double blue) {
    rasterize_rectangle(_canvas,
			x, y,
			width, height,
			Color(red, green, blue));
  }

  void vertical_semicircle(int x, int y,
			   int inner_radius, int outer_radius,
			   bool rightward,
			   double red, double green, double blue) {
    rasterize_vertical_semicircle(_canvas,
				  x, y,
				  inner_radius, outer_radius,
				  rightward,
				  Color(red, green, blue));
  }

  //////////////////////////////////////////////////////////////////////
  // Sound functions
  //////////////////////////////////////////////////////////////////////

  void make_sound(int sample_count) {
    assert(sample_count > 0);
    _sound.assign(1, sample_count, SOUND_SAMPLE_RATE, 0.0);
    assert(!is_sound_empty());
  }

  bool load_sound(const std::string& wav_file_path) {
    if (read_wav_file(_sound, wav_file_path, true)) {
      _sound.downmix();
      return true;
    } else {
      return false;
    }
  }

  bool save_sound(const std::string& wav_file_path) {
    assert(!is_sound_empty());
    return write_wav_file(_sound, wav_file_path, SOUND_BYTES_PER_SAMPLE, true);
  }

  bool is_sound_empty() {
    return _sound.is_empty();
  }

  int sound_sample_count() {
    assert(!is_sound_empty());
    return _sound.frame_count();
  }

  double get_sample(int i) {
    assert(!is_sound_empty());
    assert(i >= 0);
    assert(i < sound_sample_count());
    return _sound.get_sample(0, i);
  }

  void set_sample(int i, double sample) {
    assert(!is_sound_empty());
    assert(i >= 0);
    assert(i < sound_sample_count());
    assert(is_valid_sample(sample));
    _sound.set_sample(0, i, sample);
  }

  bool repeat_sound(double samples[],
		    int sample_count,
		    int repeat_count,
		    const std::string& wav_file_path) {
    assert(samples);
    assert(sample_count > 0);
    assert(repeat_count > 0);

    int total_samples = sample_count * repeat_count;
    Clip clip;
    clip.assign(1, total_samples, SOUND_SAMPLE_RATE, 0.0);
    
    int frame_index = 0;
    for (int i = 0; i < repeat_count; ++i)
      for (int j = 0; j < sample_count; ++j)
	clip.set_sample(0, frame_index++, samples[j]);
    assert(frame_index == total_samples);

    return write_wav_file(clip, wav_file_path, SOUND_BYTES_PER_SAMPLE, true);
  }

  //////////////////////////////////////////////////////////////////////
  // Object-oriented audio interface
  //////////////////////////////////////////////////////////////////////

  void Clip::downmix() {
    assert(!is_empty());
    if (channel_count() == 1)
      return;
    for (int i = 0; i < frame_count(); ++i) {
      double total = 0.0;
      for (int j = 0; j < channel_count(); ++j)
	total += _samples[j][i];

      double mean = total / channel_count();

      _samples[0][i] = mean;
    }

    _samples.resize(1);

    assert(!is_empty());
    assert(channel_count() == 1);
  }

  //////////////////////////////////////////////////////////////////////
  // Rasterization
  //////////////////////////////////////////////////////////////////////

  void rasterize_rectangle(Image& target,
			   int x, int y,
			   int width, int height,
			   const Color& color) {
    if ((width < 0) || (height < 0))
      return;
    for (int xi = x; xi < (x+width); ++xi)
      for (int yi = y; yi < (y+height); ++yi)
	if (target.is_valid_coordinates(xi, yi))
	  target.set_pixel(xi, yi, color);
  }

  void rasterize_vertical_semicircle(Image& target,
				     int x, int y,
				     int inner_radius, int outer_radius,
				     bool rightward,
				     const Color& color) {
    assert(inner_radius >= 0);
    assert(outer_radius >= 0);
    assert(inner_radius <= outer_radius);

    // compare the square of distances, instead of distances
    // themselves, to avoid having to call the sqrt function which is
    // expensive
    int inner2 = inner_radius*inner_radius,
        outer2 = outer_radius*outer_radius;

    // This is a simple, yet very slow, way of doing this.  It
    // suffices for now.
    for (int dx = 0; dx <= outer_radius; ++dx) {
      for (int dy = -outer_radius; dy <= outer_radius; ++dy) {
	int dist2 = dx*dx + dy*dy;
	if ((dist2 >= inner2) && (dist2 <= outer2)) {
	  int xi, yi = y + outer_radius + dy;
	  if (rightward)
	    xi = x + dx;
	  else
	    xi = x + outer_radius - dx;

	  if (target.is_valid_coordinates(xi, yi))
	    target.set_pixel(xi, yi, color);
	}
      }
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Approximate comparison
  //////////////////////////////////////////////////////////////////////

  bool approximate_equal(const Image& a, const Image& b, double delta) {
    assert(is_valid_delta(delta));

    if (a.is_empty())
      return b.is_empty();

    if ((a.get_width() != b.get_width()) || (a.get_height() != b.get_height()))
      return false;

    for (int y = 0; y < a.get_height(); ++y)
      for (int x = 0; x < a.get_width(); ++x)
	if (!approx_equal(a.get_pixel(x, y), b.get_pixel(x, y), delta))
	  return false;

    return true;
  }

  bool approximate_equal(const Clip& a, const Clip& b, double delta) {
    assert(is_valid_delta(delta));

    if (a.is_empty())
      return b.is_empty();

    if ((a.channel_count() != b.channel_count()) ||
	(a.frame_count() != b.frame_count()))
      return false;

    for (int frame = 0; frame < a.frame_count(); ++frame)
      for (int channel = 0; channel < a.channel_count(); ++channel)
	if (!approx_equal(a.get_sample(channel, frame),
			  b.get_sample(channel, frame),
			  delta))
	  return false;

    return true;
  }

  //////////////////////////////////////////////////////////////////////
  // PPM File Format
  //////////////////////////////////////////////////////////////////////

  bool write_ppm(const Image& image,
		 ostream& out,
		 bool binary,
		 bool verbose) {
    ostringstream header;
    string message;

    // file type
    message = "writing ";
    if (binary) {
      message += "binary";
      header << PPM_BINARY_HEADER;
    } else {
      message += "ASCII";
      header << PPM_ASCII_HEADER;
    }
    header << '\n';

    // dimensions
    header << image.get_width() << ' ' << image.get_height() << '\n';

    // maximum intensity value
    header << 255 << '\n';

    message += (" image, width=" + str(image.get_width()) +
		" height=" + str(image.get_height()) + " saturated=255");

    log_ppm(verbose, message);

    write_string(out, header.str());

    // pixels, row-major order
    int column = 0;
    for (int y = 0; y < image.get_height(); ++y) {
      // each row
      for (int x = 0; x < image.get_width(); ++x) {
	const Color& c = image.get_pixel(x, y);
	write_intensity(out, column, c.get_red(), binary);
	write_intensity(out, column, c.get_green(), binary);
	write_intensity(out, column, c.get_blue(), binary);
      }
    }

    // final newline to make ImageMagick happy
    out.put('\n');

    return out.good();
  }

  bool read_ppm(Image& image, istream& in, bool verbose) {
    image.clear();

    string header;
    if (!read_fixed_string(in, header, 2)) {
      log_ppm(verbose, "I/O error while reading header");
      return false;
    }

    bool binary;
    if (header == PPM_ASCII_HEADER)
      binary = false;
    else if (header == PPM_BINARY_HEADER)
      binary = true;
    else {
      log_ppm(verbose, "ppm: unrecognized file header '" + header + "'");
      return false;
    }

    int width, height, saturated;
    if (!read_int(in, width) ||
	!read_int(in, height) ||
	!read_int(in, saturated)) {
      log_ppm(verbose, "I/O error reading metadata");
      return false;
    }
    if (width <= 0) {
      log_ppm(verbose, "invalid image width " + str(width));
      return false;
    }
    if (height <= 0) {
      log_ppm(verbose, "invalid image height " + str(height));
      return false;
    }
    if (saturated <= 0) {
      log_ppm(verbose, "invalid saturation value " + str(saturated));
      return false;
    }

    log_ppm(verbose, ("found " +
		      string(binary ? "binary" : "ASCII") +
		      " image, width=" + str(width) +
		      " height=" + str(height) +
		      " saturation=" + str(saturated)));

    // one newline following the saturated digits
    char newline = in.get();
    if (!in.good()) {
      log_ppm(verbose, "I/O error immediately after header");
      return false;
    }
    if (newline != '\n') {
      log_ppm(verbose, 
	      "expected newline after metadata but found '" + newline + string("'"));
      return false;
    }

    image.resize(width, height, BLACK);

    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
	int red, green, blue;

	if (binary) {
	  red   = in.get();
	  green = in.get();
	  blue  = in.get();
	} else {
	  read_int(in, red);
	  read_int(in, green);
	  read_int(in, blue);
	}

	if (in.fail() // don't use !in.good() because an EOF is acceptable
	    || !is_in_range(0, red, saturated)
	    || !is_in_range(0, green, saturated)
	    || !is_in_range(0, blue, saturated)) {
	  log_ppm(verbose, 
		  "I/O error reading pixel at (" + str(x) + ", " + str(y) + ")");
	  image.clear();
	  return false;
	}

	image.set_pixel(x, y, Color(byte_to_intensity(red, saturated),
				    byte_to_intensity(green, saturated),
				    byte_to_intensity(blue, saturated)));
      }
    }

    return true;
  }

  bool write_ppm_file(const Image& image,
		      const std::string& path,
		      bool binary,
		      bool verbose) {
    log_ppm(verbose, "writing '" + path + "'...");
    std::ofstream f(path.c_str(), std::ofstream::binary);
    if (!f.good()) {
      log_ppm(verbose, "couldn't open file '" + path + "' for writing");
      return false;
    }
    bool ok = write_ppm(image, f, binary, verbose);
    f.close();
    return ok;
  }

  // Read a PPM file from the given path.  Inputs and outputs are the
  // same as write_ppm.
  bool read_ppm_file(Image& image,
		     const std::string& path,
		     bool verbose) {
    log_ppm(verbose, "reading '" + path + "'...");
    std::ifstream f(path.c_str(), std::ifstream::binary);
    if (!f.good()) {
      log_ppm(verbose, "couldn't open file '" + path + "' for reading");
      image.clear();
      return false;
    }
    bool ok = read_ppm(image, f, verbose);
    f.close();
    return ok;
  }

  //////////////////////////////////////////////////////////////////////
  // WAV File Format
  //////////////////////////////////////////////////////////////////////

  bool write_wav(const Clip& clip,
		 std::ostream& out,
		 int bytes_per_sample,
		 bool verbose) {
    assert(!clip.is_empty());
    assert(clip.channel_count() <= 0xFFFF); // must fit in 16 bits
    assert((bytes_per_sample == 1) || (bytes_per_sample == 2));

    // calculate metadata
    int sample_count = clip.frame_count(),
      sample_rate = clip.get_sample_rate(),
      channel_count = clip.channel_count(),
      bits_per_sample = bytes_per_sample * 8,
      bytes_per_frame = bytes_per_sample * channel_count,
      block_align = channel_count * bytes_per_sample,
      sub_chunk_2_size = sample_count * bytes_per_frame,
      chunk_size = 36 + sub_chunk_2_size;
    
    log_wav(verbose,
	    "writing channels=" + str(channel_count) +
	    ", sample-rate=" + str(sample_rate) +
	    ", bit-depth=" + str(bits_per_sample) +
	    ", samples=" + str(sample_count));

    // RIFF chunk descriptor
    if (!(write_string(out, WAV_CHUNK_ID) &&
	  write_uint(out, chunk_size, 4) &&
	  write_string(out, WAV_FORMAT))) {
      log_wav(verbose, "I/O error writing RIFF header");
      return false;
    }

    // fmt subchunk
    if (!(write_string(out, WAV_SUB_CHUNK_1_ID) &&
	  write_uint(out, 16, 4) && // subchunk size
	  write_uint(out,  1, 2) && // uncompressed linear quantization
	  write_uint(out, channel_count, 2) &&
	  write_uint(out, sample_rate, 4) &&
	  write_uint(out, sample_rate*bytes_per_frame, 4) &&
	  write_uint(out, block_align, 2) &&
	  write_uint(out, bits_per_sample, 2))) {
      log_wav(verbose, "I/O error writing fmt subchunk header");
      return false;
    }

    // data subchunk
    // two fields before the samples
    if (!(write_string(out, WAV_SUB_CHUNK_2_ID) &&
	  write_uint(out, sub_chunk_2_size, 4))) {
      log_wav(verbose, "I/O error writing data subchunk header");
      return false;
    }

    // channels are interlaced, so write each channels sample 0, then
    // each channel's sample 1, etc.
    for (int frame = 0; frame < sample_count; ++frame) {
      for (int channel = 0; channel < channel_count; ++channel) {
	double sample = clip.get_sample(channel, frame);

	unsigned bits;
	switch (bytes_per_sample) {
	case 1:
	  // [-1, +1] -> [0, 255]
	  bits = static_cast<unsigned>(((sample + 1.0) / 2.0) * 0xFF);
	  break;
	case 2:
	  // [-1, +1] -> [-32768, 32767]
	  if (sample >= 0)
	    bits = static_cast<unsigned>(sample * 0x7FFF);
	  else
	    // 2's complement
	    bits = 0xFFFF - static_cast<unsigned>(fabs(sample) * 0x8000);
	  break;
	default:
	  // shouldn't get here
	  assert(false);
	}

	if (!write_uint(out, bits, bytes_per_sample)) {
	  log_wav(verbose, "I/O error writing frame " + str(frame));
	  return false;
	}
      }
    }

    return true;
  }

  bool read_wav(Clip& clip,
		std::istream& in,
		bool verbose) {
    clip.clear();

    // RIFF header
    string id, format;
    unsigned chunk_size;
    if (!(read_fixed_string(in, id, 4) &&
	  read_uint(in, chunk_size, 4) &&
	  read_fixed_string(in, format, 4))) {
      log_wav(verbose, "I/O error reading RIFF header");
      return false;
    }
    if (id != WAV_CHUNK_ID) {
      log_wav(verbose, "invalid chunk ID '" + id + "'");
      return false;
    }
    if (format != WAV_FORMAT) {
      log_wav(verbose, "invalid format subchunk id '" + id + "'");
      return false;
    }

    // fmt subchunk
    unsigned subchunk_1_size,
      audio_format,
      num_channels,
      sample_rate,
      byte_rate, 
      block_align,
      bits_per_sample;
    if (!(read_fixed_string(in, id, 4) &&
	  read_uint(in, subchunk_1_size, 4) &&
	  read_uint(in, audio_format, 2) &&
	  read_uint(in, num_channels, 2) &&
	  read_uint(in, sample_rate, 4) &&
	  read_uint(in, byte_rate, 4) &&
	  read_uint(in, block_align, 2) &&
	  read_uint(in, bits_per_sample, 2))) {
      log_wav(verbose, "I/O error in format subchunk 1 header");
      return false;
    }
    if (id != WAV_SUB_CHUNK_1_ID) {
      log_wav(verbose, "invalid subchunk 1 id '" + id + "'");
      return false;
    }
    if ((subchunk_1_size != 16) && verbose)
      log_wav(verbose, "ignoring malformed subchunk 1 size " + str(subchunk_1_size));
    if (audio_format != 1) {
      log_wav(verbose, "unsupported audio format " + str(audio_format));
      return false;
    }
    if (num_channels == 0) {
      log_wav(verbose, "invalid channel count zero");
      return false;
    }
    if (sample_rate == 0) {
      log_wav(verbose, "invalid sample rate zero");
      return false;
    }
    // check bits per sample first since the next two depend on it
    if ((bits_per_sample != 8) && (bits_per_sample != 16)) {
      log_wav(verbose, "unsupported bit depth " + str(bits_per_sample));
      return false;
    }
    // back in order
    assert((bits_per_sample % 8) == 0);
    unsigned bytes_per_sample = bits_per_sample / 8;
    unsigned expected_byte_rate = sample_rate * num_channels * bytes_per_sample;
    if (byte_rate != expected_byte_rate) {
      log_wav(verbose, "correcting malformed byte rate " + str(byte_rate) + " to " + str(expected_byte_rate));
      byte_rate = expected_byte_rate;
    }
    unsigned expected_block_align = num_channels * bytes_per_sample;
    if (block_align != expected_block_align) {
      log_wav(verbose, "correcting malformed block-align " + str(block_align) + " to " + str(expected_block_align));
      block_align = expected_block_align;
    }

    // data subchunk
    unsigned subchunk_2_size;
    if (!(read_fixed_string(in, id, 4) &&
	  read_uint(in, subchunk_2_size, 4))) {
      log_wav(verbose, "I/O error reading data subchunk header");
      return false;
    }
    if (id != WAV_SUB_CHUNK_2_ID) {
      log_wav(verbose, "invalid subchunk 2 id '" + id + "'");
      return false;
    }

    unsigned num_samples = subchunk_2_size / (num_channels * bytes_per_sample);

    log_wav(verbose,
	    "found channels=" + str(num_channels) +
	    ", sample-rate=" + str(sample_rate) +
	    ", bit-depth=" + str(bits_per_sample) +
	    ", samples=" + str(num_samples));

    clip.assign(num_channels, num_samples, sample_rate, 0.0);

    // each frame
    for (unsigned frame = 0; frame < num_samples; ++frame) {
      // each channel
      for (unsigned channel = 0; channel < num_channels; ++channel) {
	unsigned bits;
	if (!read_uint(in, bits, bytes_per_sample)) {
	  log_wav(verbose, "I/O error writing frame " + str(frame));
	  clip.clear();
	  return false;
	}

	double sample;
	switch (bytes_per_sample) {
	case 1:
	  sample = ((bits / 255.0) * 2.0) - 1.0;
	  break;
	case 2:
	  if ((bits & (1 << 15)) == 0) {
	    // positive
	    sample = bits / static_cast<double>(0x7FFF);
	    assert(is_valid_sample(sample));
	  } else {
	    // twos complement negative
	    sample = -static_cast<double>(0xFFFF - bits) / 0x8000;
	    assert(is_valid_sample(sample));
	  }
	  break;
	default:
	  assert(false);
	}

	clip.set_sample(channel, frame, sample);
      }
    }

    return true;
  }

  bool write_wav_file(const Clip& clip,
		      const std::string& path,
		      int bytes_per_sample,
		      bool verbose) {
    log_wav(verbose, "writing '" + path + "'...");
    std::ofstream f(path.c_str(), std::ofstream::binary);
    if (!f.good()) {
      log_wav(verbose, "couldn't open file '" + path + "' for writing");
      return false;
    }
    bool ok = write_wav(clip, f, bytes_per_sample, verbose);
    f.close();
    return ok;
  }

  bool read_wav_file(Clip& clip,
		     const std::string& path,
		     bool verbose) {
    log_wav(verbose, "reading '" + path + "'...");
    std::ifstream f(path.c_str(), std::ifstream::binary);
    if (!f.good()) {
      log_wav(verbose, "couldn't open file '" + path + "' for reading");
      clip.clear();
      return false;
    }
    bool ok = read_wav(clip, f, verbose);
    f.close();
    return ok;
  }

} // namespace csuf
