
//////////////////////////////////////////////////////////////////////
// csufmedia_test.cpp
//
// See csufmedia.h for substantive comments.
//
// See LICENSE for copyright and license information.
// 
//////////////////////////////////////////////////////////////////////

#include <cmath>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "csufmedia.h"

using namespace csuf;
using namespace std;

const double EPS = 1e-3,  // tolerance for floating point comparison
             DELTA = .10;
const string PPM_FILENAME = "test.ppm",
             PNG_FILENAME = "test.png";

const int FRAME_COUNT = 1000000,
  SAMPLE_RATE = 44100,
  NOTE_FREQ = 440;
const string WAV_FILENAME = "test.wav",
  WAV_2_FILENAME = "test2.wav",
  FLAC_FILENAME = "test.flac";

BOOST_AUTO_TEST_CASE( test_is_in_range ) {
  // int overload
  BOOST_CHECK(is_in_range(1, 5, 10));   // yes
  BOOST_CHECK(!is_in_range(1, -5, 10)); // no
  BOOST_CHECK(is_in_range(1, 1, 10));   // inclusive
  BOOST_CHECK(is_in_range(1, 10, 10));  // inclusive

  // double overload
  BOOST_CHECK(is_in_range(1.0, 5.0, 10.0));   // yes
  BOOST_CHECK(!is_in_range(1.0, -5.0, 10.0)); // no
  BOOST_CHECK(is_in_range(1.0, 1.0, 10.0));   // inclusive
  BOOST_CHECK(is_in_range(1.0, 10.0, 10.0));  // inclusive
}

BOOST_AUTO_TEST_CASE( test_is_valid_intensity ) {
  BOOST_CHECK(is_valid_intensity(0.0));
  BOOST_CHECK(is_valid_intensity(0.5));
  BOOST_CHECK(is_valid_intensity(1.0));
  BOOST_CHECK(!is_valid_intensity(-1.0));
  BOOST_CHECK(!is_valid_intensity(1.01));
}

BOOST_AUTO_TEST_CASE( test_Color ) {
  // constructor and getters
  Color red(1, 0, 0);
  BOOST_CHECK_EQUAL(1.0, red.get_red());
  BOOST_CHECK_EQUAL(0.0, red.get_green());
  BOOST_CHECK_EQUAL(0.0, red.get_blue());
  Color green(0, 1, 0);
  BOOST_CHECK_EQUAL(0.0, green.get_red());
  BOOST_CHECK_EQUAL(1.0, green.get_green());
  BOOST_CHECK_EQUAL(0.0, green.get_blue());
  Color blue(0, 0, 1);
  BOOST_CHECK_EQUAL(0.0, blue.get_red());
  BOOST_CHECK_EQUAL(0.0, blue.get_green());
  BOOST_CHECK_EQUAL(1.0, blue.get_blue());

  // Mutators
  const double before = 0, after = 0.6;
  Color c = Color(before, before, before);
  // set_red
  BOOST_CHECK_CLOSE(before, c.get_red(), EPS);
  c.set_red(after);
  BOOST_CHECK_CLOSE(after, c.get_red(), EPS);
  // set_green
  BOOST_CHECK_CLOSE(before, c.get_green(), EPS);
  c.set_green(after);
  BOOST_CHECK_CLOSE(after, c.get_green(), EPS);
  // set_blue
  BOOST_CHECK_CLOSE(before, c.get_blue(), EPS);
  c.set_blue(after);
  BOOST_CHECK_CLOSE(after, c.get_blue(), EPS);
}

BOOST_AUTO_TEST_CASE( test_color_constants ) {
  BOOST_CHECK_EQUAL(0, BLACK.get_red());
  BOOST_CHECK_EQUAL(0, BLACK.get_green());
  BOOST_CHECK_EQUAL(0, BLACK.get_blue());

  BOOST_CHECK_EQUAL(1, WHITE.get_red());
  BOOST_CHECK_EQUAL(1, WHITE.get_green());
  BOOST_CHECK_EQUAL(1, WHITE.get_blue());

  BOOST_CHECK_EQUAL(1, RED.get_red());
  BOOST_CHECK_EQUAL(0, RED.get_green());
  BOOST_CHECK_EQUAL(0, RED.get_blue());

  BOOST_CHECK_EQUAL(0, GREEN.get_red());
  BOOST_CHECK_EQUAL(1, GREEN.get_green());
  BOOST_CHECK_EQUAL(0, GREEN.get_blue());

  BOOST_CHECK_EQUAL(0, BLUE.get_red());
  BOOST_CHECK_EQUAL(0, BLUE.get_green());
  BOOST_CHECK_EQUAL(1, BLUE.get_blue());
}

BOOST_AUTO_TEST_CASE( test_Image_empty ) {
  // constructor
  Image im;
  BOOST_CHECK(im.is_empty());
  const int BIG = 200;
  for (int x = -BIG; x <= BIG; ++x)
    for (int y = -BIG; y <= BIG; ++y)
      BOOST_CHECK(!im.is_valid_coordinates(x, y));

  // operator==
  Image im2;
  BOOST_CHECK(im == im2);

  // clear
  Image clear;
  BOOST_CHECK(clear.is_empty());
  im.clear();
  BOOST_CHECK(clear.is_empty());

  // resize
  Image resize;
  BOOST_CHECK(resize.is_empty());
  resize.resize(100, 200, BLACK);
  BOOST_CHECK(!resize.is_empty());
  BOOST_CHECK_EQUAL(100, resize.get_width());
  BOOST_CHECK_EQUAL(200, resize.get_height());  

  // fill
  Image fill;
  BOOST_CHECK(fill.is_empty());
  fill.fill(BLACK);
  BOOST_CHECK(fill.is_empty());
}

BOOST_AUTO_TEST_CASE( test_Image_nonempty ) {
  // sized constructor
  // is_empty
  // get_width
  // get_height
  const int W = 100, H = 200;
  BOOST_CHECK(W != H); // should differ
  Image im(W, H, BLACK);
  BOOST_CHECK(!im.is_empty());
  BOOST_CHECK_EQUAL(W, im.get_width());  // get_width
  BOOST_CHECK_EQUAL(H, im.get_height()); // get_height
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y)
      BOOST_CHECK(im.get_pixel(x, y) == BLACK);

  // operator==
  // identical
  Image im2(W, H, BLACK);
  BOOST_CHECK(im == im2);
  // same dimensions, one different pixel
  Image im3 = im2;
  im3.set_pixel(W/2, H/2, RED);
  BOOST_CHECK(!(im2 == im3));
  // same colors, different dimensions
  Image im4(H, W, BLACK);
  BOOST_CHECK(!(im3 == im4));

  // is_valid_coordinates
  BOOST_CHECK(im.is_valid_coordinates(0, 0));
  BOOST_CHECK(im.is_valid_coordinates(W-1, H-1));
  BOOST_CHECK(!im.is_valid_coordinates(-1, -1));
  BOOST_CHECK(!im.is_valid_coordinates(W, H));
  BOOST_CHECK(!im.is_valid_coordinates(W*2, H*2));
  const int BIG = 2 * max(W, H);
  for (int x = -BIG; x <= BIG; ++x)
    for (int y = -BIG; y <= BIG; ++y)
      BOOST_CHECK_EQUAL((x >= 0) && (y >= 0) && (x < W) && (y < H),
			im.is_valid_coordinates(x, y) );

  // clear
  Image cleared;
  BOOST_CHECK(cleared.is_empty());
  cleared.resize(W, H, BLACK);
  BOOST_CHECK(!cleared.is_empty());

  // resize
  Image resized;
  BOOST_CHECK(resized.is_empty());
  // from empty
  resized.resize(10, 11, BLACK);
  BOOST_CHECK(!resized.is_empty());
  BOOST_CHECK_EQUAL(10, resized.get_width()); 
  BOOST_CHECK_EQUAL(11, resized.get_height()); 
  BOOST_CHECK(resized.get_pixel(0, 0) == BLACK);
  // enlarge
  resized.resize(20, 21, RED);
  BOOST_CHECK(!resized.is_empty());
  BOOST_CHECK_EQUAL(20, resized.get_width()); 
  BOOST_CHECK_EQUAL(21, resized.get_height()); 
  BOOST_CHECK(resized.get_pixel(0, 0) == RED);

  // fill
  Image fill(W, H, BLACK);
  fill.fill(RED);
  BOOST_CHECK(fill.get_pixel(0, 0) == RED);
  BOOST_CHECK(fill.get_pixel(W-1, H-1) == RED);
  fill.fill(WHITE);
  BOOST_CHECK(fill.get_pixel(0, 0) == WHITE);
  BOOST_CHECK(fill.get_pixel(W-1, H-1) == WHITE);

  // get_pixel
  // get_mutable_pixel
  Image mutate(W, H, BLACK);
  BOOST_CHECK(mutate.get_pixel(0, 0) == BLACK);
  Color& pixel = mutate.get_mutable_pixel(0, 0);
  pixel.set_red(1);
  BOOST_CHECK(mutate.get_pixel(0, 0) == RED);

  // set_pixel
  Image set(W, H, BLACK);
  BOOST_CHECK(set.get_pixel(0, 0) == BLACK);
  set.set_pixel(0, 0, GREEN);
  BOOST_CHECK(set.get_pixel(0, 0) == GREEN);
  
}

BOOST_AUTO_TEST_CASE( test_read_ppm_file ) {
  // nonexistent file
  Image junk;
  BOOST_CHECK(!read_ppm_file(junk, "<NOTFOUND>.ppt", false));
  BOOST_CHECK(junk.is_empty());

  // malformed file
  BOOST_CHECK(!read_ppm_file(junk, "SConstruct", false));
  BOOST_CHECK(junk.is_empty());

  Image gimp_ascii;
  BOOST_CHECK(read_ppm_file(gimp_ascii, "elephant_gimp_ascii.ppm", true));
  BOOST_CHECK(!gimp_ascii.is_empty());
  Image proper = gimp_ascii;
  BOOST_CHECK(gimp_ascii == proper);
  BOOST_CHECK(approximate_equal(gimp_ascii, proper, DELTA));

  Image gimp_binary;
  BOOST_CHECK(read_ppm_file(gimp_binary, "elephant_gimp_binary.ppm", true));
  BOOST_CHECK(approximate_equal(gimp_binary, proper, DELTA));

  Image irfan_ascii;
  BOOST_CHECK(read_ppm_file(irfan_ascii, "elephant_irfan_ascii.ppm", true));
  BOOST_CHECK(approximate_equal(irfan_ascii, proper, DELTA));

  Image irfan_binary;
  BOOST_CHECK(read_ppm_file(irfan_binary, "elephant_irfan_ascii.ppm", true));
  BOOST_CHECK(approximate_equal(irfan_binary, proper, DELTA));

  Image magick;
  BOOST_CHECK(read_ppm_file(magick, "elephant_magick.ppm", true));
  BOOST_CHECK(approximate_equal(magick, proper, DELTA));
}

BOOST_AUTO_TEST_CASE( test_write_ppm_ascii ) {
  Image orig;
  BOOST_CHECK(read_ppm_file(orig, "elephant_gimp_ascii.ppm", true));
  BOOST_CHECK(write_ppm_file(orig, PPM_FILENAME, false, true));
  Image fresh;
  BOOST_CHECK(read_ppm_file(fresh, PPM_FILENAME, true));
  BOOST_CHECK(fresh == orig);
  remove(PPM_FILENAME.c_str());
}

BOOST_AUTO_TEST_CASE( test_write_ppm_binary ) {
  Image orig;
  BOOST_CHECK(read_ppm_file(orig, "elephant_gimp_binary.ppm", true));
  BOOST_CHECK(write_ppm_file(orig, PPM_FILENAME, true, true));
  Image fresh;
  BOOST_CHECK(read_ppm_file(fresh, PPM_FILENAME, true));
  BOOST_CHECK(fresh == orig);
  remove(PPM_FILENAME.c_str());
}

BOOST_AUTO_TEST_CASE( test_ppm_ascii_round_trip_magick_compatible ) {
  // write a fresh ASCII PPM
  Image orig;
  BOOST_CHECK(read_ppm_file(orig, "elephant_gimp_binary.ppm", true));
  BOOST_CHECK(write_ppm_file(orig, PPM_FILENAME, false, true));

  // have ImageMagick convert our PPM to PNG
  string command = "convert " + PPM_FILENAME + " " + PNG_FILENAME;
  BOOST_CHECK_EQUAL(0, system(command.c_str()));

  // convert its PNG to PPM
  command = "convert " + PNG_FILENAME + " " + PPM_FILENAME;
  BOOST_CHECK_EQUAL(0, system(command.c_str()));

  // load ImageMagick's PPM
  Image result;
  BOOST_CHECK(read_ppm_file(result, PPM_FILENAME, true));

  // images identical?
  BOOST_CHECK(result == orig);

  // delete temporary files
  remove(PNG_FILENAME.c_str());
  remove(PPM_FILENAME.c_str());
}

BOOST_AUTO_TEST_CASE( test_ppm_binary_round_trip_magick_compatible ) {
  Image orig;
  BOOST_CHECK(read_ppm_file(orig, "elephant_gimp_binary.ppm", true));
  BOOST_CHECK(write_ppm_file(orig, PPM_FILENAME, true, true));

  string command = "convert " + PPM_FILENAME + " " + PNG_FILENAME;
  BOOST_CHECK_EQUAL(0, system(command.c_str()));

  command = "convert " + PNG_FILENAME + " " + PPM_FILENAME;
  BOOST_CHECK_EQUAL(0, system(command.c_str()));

  Image result;
  BOOST_CHECK(read_ppm_file(result, PPM_FILENAME, true));

  BOOST_CHECK(result == orig);

  remove(PNG_FILENAME.c_str());
  remove(PPM_FILENAME.c_str());
}

inline void sin_clip(Clip& clip, int channel_count) {
  double period = SAMPLE_RATE / NOTE_FREQ;
  clip.assign(channel_count, FRAME_COUNT, 44100, 0.0);
  for (int frame = 0; frame < FRAME_COUNT; ++frame) {
    for (int channel = 0; channel < channel_count; ++channel) {
      double radians = (2.0*M_PI*frame/period) + (channel*M_PI/4.0);
      clip.set_sample(channel, frame, .9 * sin(radians));
    }
  }
}

BOOST_AUTO_TEST_CASE( test_read_write_wav ) {
  for (int channels = 1; channels <= 2; ++channels) {
    for (int bytes_per_sample = 1; bytes_per_sample <= 2; ++bytes_per_sample) {
      Clip before;

      sin_clip(before, channels);

      BOOST_ASSERT(write_wav_file(before, WAV_FILENAME, bytes_per_sample, true));

      Clip after;
      BOOST_ASSERT(read_wav_file(after, WAV_FILENAME, true));

      remove(WAV_FILENAME.c_str());

      BOOST_ASSERT(approximate_equal(before, after, DELTA));
    }
  }
}

BOOST_AUTO_TEST_CASE( test_wav_sox_round_trip ) {
  for (int channels = 1; channels <= 2; ++channels) {
    for (int bytes_per_sample = 1; bytes_per_sample <= 2; ++bytes_per_sample) {
      Clip before;

      sin_clip(before, channels);

      BOOST_ASSERT(write_wav_file(before, WAV_FILENAME, bytes_per_sample, true));

      string command = "sox " + WAV_FILENAME + " " + FLAC_FILENAME;
      system(command.c_str());
      command = "sox " + FLAC_FILENAME + " " + WAV_2_FILENAME;
      system(command.c_str());

      Clip after;
      BOOST_ASSERT(read_wav_file(after, WAV_2_FILENAME, true));

      remove(WAV_FILENAME.c_str());
      remove(WAV_2_FILENAME.c_str());
      remove(FLAC_FILENAME.c_str());

      BOOST_ASSERT(approximate_equal(before, after, DELTA));
    }
  }
}
