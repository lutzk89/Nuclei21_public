/*
LodePNG Examples

Copyright (c) 2005-2010 Lode Vandevenne

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software
    in a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

    2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

    3. This notice may not be removed or altered from any source
    distribution.
*/

/*These examples demonstrate the usage of LodePNG in C++ and are useful as test.
Place the example code in a new .cpp file, and compile together with lodepng.cpp
and lodepng.h, with optimization.
Run the example executables from command line with the filename as parameter.*/

//g++ *.cpp -lSDL -Wall -Wextra -pedantic -ansi -lSDL -O3


////////////////////////////////////////////////////////////////////////////////
// LodePNG console example
// This example decodes a PNG file and displays its size in the console
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include "lodepng.h"

void example1(int argc, char *argv[])
{
  const char* filename = argc > 1 ? argv[1] : "test.png";
  
  //load and decode
  std::vector<unsigned char> buffer, image;
  LodePNG::loadFile(buffer, filename); //load the image file with given filename
  LodePNG::Decoder decoder;
  decoder.decode(image, buffer.empty() ? 0 : &buffer[0], (unsigned)buffer.size()); //decode the png
  
  //if there's an error, display it, otherwise display information about the image
  if(decoder.hasError()) std::cout << "error " << decoder.getError() << ": " << LodePNG_error_text(decoder.getError()) << std::endl;
  else
  {
    std::cout << "\n" <<
      "w: " << decoder.getWidth() << "\n" <<
      "h: " << decoder.getHeight() << "\n" <<
      "bitDepth: " << decoder.getInfoPng().color.bitDepth << "\n" <<
      "bpp: " << decoder.getBpp() << "\n" <<
      "colorChannels: " << decoder.getChannels() << "\n" <<
      "paletteSize: " << decoder.getInfoPng().color.palettesize / 4 << "\n" <<
      "colorType: " << decoder.getInfoPng().color.colorType << "\n" <<
      "compressionMethod: " << decoder.getInfoPng().compressionMethod << "\n" <<
      "filterMethod: " << decoder.getInfoPng().filterMethod << "\n" <<
      "interlaceMethod: " << decoder.getInfoPng().interlaceMethod << "\n";
    for(size_t i = 0; i < decoder.getInfoPng().text.num; i++)
      std::cout << decoder.getInfoPng().text.keys[i] << ": " << decoder.getInfoPng().text.strings[i] << "\n";
    if(decoder.infoPng.time_defined)
    {
      std::cout << "modification time: " << decoder.infoPng.time.year << "-"
                << std::setw(2) << std::setfill('0') << (int)decoder.infoPng.time.month << "-"
                << std::setw(2) << std::setfill('0') << (int)decoder.infoPng.time.day << " "
                << std::setw(2) << std::setfill('0') << (int)decoder.infoPng.time.hour << ":"
                << std::setw(2) << std::setfill('0') << (int)decoder.infoPng.time.minute << ":"
                << std::setw(2) << std::setfill('0') << (int)decoder.infoPng.time.second << std::endl;
    }
    if(decoder.infoPng.phys_defined)
    {
      std::cout << "physical size: " << decoder.infoPng.phys_x << " " << decoder.infoPng.phys_y << " "
                << (int)decoder.infoPng.phys_unit << std::endl;
    }
  }
}

#define ENABLE_SDL_EXAMPLE
#ifdef ENABLE_SDL_EXAMPLE

////////////////////////////////////////////////////////////////////////////////
// LodePNG SDL example
// This example displays a PNG with a checkerboard pattern to show tranparency
// It requires the SDL library to compile and run
// g++ *.cpp -lSDL -O3
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <SDL/SDL.h> //requires SDL
#include "lodepng.h"

void example2(int argc, char *argv[])
{
  const char* filename = argc > 1 ? argv[1] : "test.png";
  
  std::vector<unsigned char> buffer, image;
  LodePNG::loadFile(buffer, filename); //load the image file with given filename
  LodePNG::Decoder decoder;
  decoder.decode(image, buffer); //decode the png
  
  //get width, height and pixels
  unsigned w = decoder.getWidth(), h =  decoder.getHeight();

  //stop if there is an error
  if(decoder.hasError())
  {
    std::cout << "error " << decoder.getError() << ": " << LodePNG_error_text(decoder.getError()) << std::endl;
    return;
  }
  
  //avoid too large window size by downscaling large image
  unsigned jump = 1;
  if(w / 1024 >= jump) jump = w / 1024 + 1;
  if(h / 1024 >= jump) jump = h / 1024 + 1;
  
  //init SDL
  if(SDL_Init(SDL_INIT_VIDEO) < 0) return;
  SDL_Surface* scr = SDL_SetVideoMode(w / jump, h / jump, 32, SDL_HWSURFACE);
  if(!scr) return;
  SDL_WM_SetCaption(filename, NULL); //set window caption
  
  //plot the pixels of the PNG file
  for(unsigned y = 0; y + jump - 1 < h; y += jump)
  for(unsigned x = 0; x + jump - 1 < w; x += jump)
  {
    //get RGBA components
    Uint32 r = image[4 * y * w + 4 * x + 0]; //red
    Uint32 g = image[4 * y * w + 4 * x + 1]; //green
    Uint32 b = image[4 * y * w + 4 * x + 2]; //blue
    Uint32 a = image[4 * y * w + 4 * x + 3]; //alpha
    
    //make translucency visible by placing checkerboard pattern behind image
    int checkerColor = 191 + 64 * (((x / 16) % 2) == ((y / 16) % 2));
    r = (a * r + (255 - a) * checkerColor) / 255;
    g = (a * g + (255 - a) * checkerColor) / 255;
    b = (a * b + (255 - a) * checkerColor) / 255;
    
    //give the color value to the pixel of the screenbuffer
    Uint32* bufp;
    bufp = (Uint32 *)scr->pixels + (y * scr->pitch / 4) / jump + (x / jump);
    *bufp = 65536 * r + 256 * g + b;
  }
  
  //pause until you press escape and meanwhile redraw screen
  SDL_Event event;
  int done = 0;
  while(done == 0)
  {
    while(SDL_PollEvent(&event))
    {
      if(event.type == SDL_QUIT) done = 1;
      if(SDL_GetKeyState(NULL)[SDLK_ESCAPE]) done = 1;
    }
    SDL_UpdateRect(scr, 0, 0, 0, 0); //redraw screen
    SDL_Delay(5); //pause 5 ms so it consumes less processing power
  }
  
  SDL_Quit();
}

#endif


////////////////////////////////////////////////////////////////////////////////
// LodePNG Encoder example
// This example encodes a PNG containing a generated texture
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "lodepng.h"

void example3(int argc, char *argv[])
{
  //check if user gave a filename
  if(argc < 2)
  {
    std::cout << "please provide a filename to save to" << std::endl;
    return;
  }

  //generate some image
  const int w = 512;
  const int h = 512;
  std::vector<unsigned char> image;
  image.resize(w * h * 4);
  for(int y = 0; y < h; y++)
  for(int x = 0; x < w; x++)
  {
    //pattern 1
    image[4 * w * y + 4 * x + 0] = (unsigned char)(127 * (1 + std::sin((                    x * x +                     y * y) / (w * h / 8.0))));
    image[4 * w * y + 4 * x + 1] = (unsigned char)(127 * (1 + std::sin(((w - x - 1) * (w - x - 1) +                     y * y) / (w * h / 8.0))));
    image[4 * w * y + 4 * x + 2] = (unsigned char)(127 * (1 + std::sin((                    x * x + (h - y - 1) * (h - y - 1)) / (w * h / 8.0))));
    image[4 * w * y + 4 * x + 3] = (unsigned char)(127 * (1 + std::sin(((w - x - 1) * (w - x - 1) + (h - y - 1) * (h - y - 1)) / (w * h / 8.0))));
    
    //pattern 2
    //image[4 * w * y + 4 * x + 0] = 255 * !(x & y);
    //image[4 * w * y + 4 * x + 1] = x ^ y;
    //image[4 * w * y + 4 * x + 2] = x | y;
    //image[4 * w * y + 4 * x + 3] = 255;
  }
  
  //create encoder and set settings and info (optional)
  LodePNG::Encoder encoder;
  encoder.addText("Comment", "Created with LodePNG");
  encoder.getSettings().zlibsettings.windowSize = 2048;

  //encode and save
  std::vector<unsigned char> buffer;
  encoder.encode(buffer, image.empty() ? 0 : &image[0], w, h);
  LodePNG::saveFile(buffer, argv[1]);
}


////////////////////////////////////////////////////////////////////////////////
// LodePNG greyscale example
// Load a PNG directly as greyscale if it's greyscale, convert otherwise
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "lodepng.h"

void example4(int argc, char *argv[])
{
  const char* filename = argc > 1 ? argv[1] : "test.png";
  
  //load and decode
  std::vector<unsigned char> buffer; //the file
  LodePNG::loadFile(buffer, filename); //load the image file with given filename
  LodePNG::Decoder decoder;
  
  std::vector<unsigned char> grey; //in here the greyscale image will be stored
  
  decoder.inspect(buffer); //get information from the header, to see if it's greyscale
  
  if(!decoder.hasError())
  {
    if(decoder.isGreyscaleType()) //only if the PNG image is greyscale, can the output be converted to greyscale by LodePNG
    {
      std::cout << "it's a greyscale PNG\n";
      decoder.getInfoRaw().color.colorType = 0; //set color type "greyscale" for the output
      decoder.decode(grey, buffer);
    }
    else //else you have to do color math to convert to greyscale
    {
      std::cout << "it's a color PNG\n";
      decoder.getInfoRaw().color.colorType = 2; //set color type "RGB" for the output
      std::vector<unsigned char> rgb; //temporary rgb image
      decoder.decode(rgb, buffer);
      
      grey.resize(rgb.size() / 3);
      for(size_t i = 0; i < grey.size(); i++) grey[i] = (rgb[i * 3 + 0] + rgb[i * 3 + 1] + rgb[i * 3 + 2]) / 3; //rgb to greyscale by taking average
    }
  }
  
  if(decoder.hasError()) std::cout << "error " << decoder.getError() << ": " << LodePNG_error_text(decoder.getError()) << std::endl;
  else std::cout << "read " << grey.size() << " greyscale pixels\n";
}


////////////////////////////////////////////////////////////////////////////////
// LodePNG palette example
// This example encodes a PNG with a palette
// Both the image and palette contain sine waves, resulting in a sort of plasma
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "lodepng.h"

void example5(int argc, char *argv[])
{
  //check if user gave a filename
  if(argc < 2)
  {
    std::cout << "please provide a filename to save to" << std::endl;
    return;
  }

  //create encoder and set settings and info (optional)
  LodePNG::Encoder encoder;
  encoder.addText("Comment", "Created with LodePNG");
  
  //Generate palette. Note how only the palette of infoPng has to be set, not that of infoRaw. That is because
  //we'll encode the PNG as a palettized PNG and the palette stored in InfoPng is used for the PLTE and tRNS
  //chunk of the PNG. When encoding a palettized raw image to a PNG without palette (RGB, RGBA, ...), you'd
  //have to set the palette of infoRaw instead of that of infoPng, because then the palette info of infoRaw is
  //used by the encoder to do the conversion to RGBA.
  for(int i = 0; i < 256; i++)
  {
    encoder.addPalette((unsigned char)(127 * (1 + sin(5 * i * 6.28318531 / 256)))
                     , (unsigned char)(127 * (1 + sin(2 * i * 6.28318531 / 256)))
                     , (unsigned char)(127 * (1 + sin(3 * i * 6.28318531 / 256)))
                     , (unsigned char)( 63 * (1 + sin(8 * i * 6.28318531 / 256))) + 128); /*alpha channel of the palette (tRNS chunk)*/
  }
  
  //both the raw image and the encoded image must get colorType 3 (palette)
  encoder.getInfoPng().color.colorType = 3; //if you comment this line, and store the palette in InfoRaw instead (use getInfoRaw() in the previous lines), then you get the same image in a RGBA PNG.
  encoder.getInfoRaw().color.colorType = 3;
  
  //generate some image
  const unsigned w = 512;
  const unsigned h = 512;
  std::vector<unsigned char> image;
  image.resize(w * h);
  for(unsigned y = 0; y < h; y++)
  for(unsigned x = 0; x < w; x++)
  {
    image[w * y + x] = (unsigned char)(63 * ((1 + std::sin((2 * 6.28318531) * (x / float(w)))) + (1 + std::sin((2 * 6.28318531) * (y / float(h))))));
  }
  
  //encode and save
  std::vector<unsigned char> buffer;
  encoder.encode(buffer, image.empty() ? 0 : &image[0], w, h);
  if(encoder.hasError()) std::cout << "Encoder error " << encoder.getError() << ": " << LodePNG_error_text(encoder.getError()) << std::endl;
  LodePNG::saveFile(buffer, argv[1]);
}

////////////////////////////////////////////////////////////////////////////////
// LodePNG Ghostifier
// Loads a PNG, modifies its alpha channel to a ghost effect, saves the result
// The result is visible in browsers that support PNG transparency
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "lodepng.h"

void example6(int argc, char *argv[])
{
  //check if user gave a filename
  if(argc < 2)
  {
    std::cout << "please provide a filename to save to" << std::endl;
    return;
  }
  std::string filename = argv[1];
  
  std::vector<unsigned char> image;
  unsigned w, h;
  
  int error = LodePNG::decode(image, w, h, filename);

  if(error) return;
  
  //ghostify image
  for(size_t y = 0; y < h; y++)
  for(size_t x = 0; x < w; x++)
  {
    size_t index = 4 * w * y + 4 * x;
    image[index + 3] = (unsigned char)((int(image[index + 0]) + int(image[index + 1]) + int(image[index + 2])) / 3);
  }
  
  //encode and save
  LodePNG::encode("ghostified_" + filename, image.empty() ? 0 : &image[0], w, h);
}

////////////////////////////////////////////////////////////////////////////////
// LodePNG 4-bit palette example
// This example encodes a 511x511 PNG with a 4-bit palette
// Both image and palette contain sine waves, resulting in a sort of plasma
////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "lodepng.h"

void example7(int argc, char *argv[])
{
  //check if user gave a filename
  if(argc < 2)
  {
    std::cout << "please provide a filename to save to" << std::endl;
    return;
  }

  //create encoder and set settings and info (optional)
  LodePNG::Encoder encoder;
  encoder.addText("Comment", "Created with LodePNG");
  
  //Generate palette. Note how only the palette of infoPng has to be set, not that of infoRaw. That is because
  //we'll encode the PNG as a palettized PNG and the palette stored in InfoPng is used for the PLTE and tRNS
  //chunk of the PNG. When encoding a palettized raw image to a PNG without palette (RGB, RGBA, ...), you'd
  //have to set the palette of infoRaw instead of that of infoPng, because then the palette info of infoRaw is
  //used by the encoder to do the conversion to RGBA.
  for(int i = 0; i < 16; i++)
  {
    encoder.addPalette((unsigned char)(127 * (1 + std::sin(5 * i * 6.28318531 / 16)))
                     , (unsigned char)(127 * (1 + std::sin(2 * i * 6.28318531 / 16)))
                     , (unsigned char)(127 * (1 + std::sin(3 * i * 6.28318531 / 16)))
                     , (unsigned char)( 63 * (1 + std::sin(8 * i * 6.28318531 / 16))) + 128); /*alpha channel of the palette (tRNS chunk)*/
  }
  
  //both the raw image and the encoded image must get colorType 3 (palette)
  encoder.getInfoPng().color.colorType = 3; //if you comment this line, and store the palette in InfoRaw instead (use getInfoRaw() in the previous lines), then you get the same image in a RGBA PNG.
  encoder.getInfoRaw().color.colorType = 3;
  encoder.getInfoPng().color.bitDepth = 4;
  encoder.getInfoRaw().color.bitDepth = 4;
  
  //generate some image
  const unsigned w = 511;
  const unsigned h = 511;
  std::vector<unsigned char> image;
  image.resize((w * h * 4 + 7) / 8, 0);
  for(unsigned y = 0; y < h; y++)
  for(unsigned x = 0; x < w; x++)
  {
    size_t byte_index = (y * w + x) / 2;
    bool byte_half = (y * w + x) % 2 == 1;
    
    int color = (int)(4 * ((1 + std::sin(2.0 * 6.28318531 * x / (double)w))
                         + (1 + std::sin(2.0 * 6.28318531 * y / (double)h))) );

    image[byte_index] |= (unsigned char)(color << (byte_half ? 0 : 4));
  }
  
  //encode and save
  std::vector<unsigned char> buffer;
  encoder.encode(buffer, image.empty() ? 0 : &image[0], w, h);
  if(encoder.hasError()) std::cout << "Encoder error " << encoder.getError() << ": " << LodePNG_error_text(encoder.getError()) << std::endl;
  LodePNG::saveFile(buffer, argv[1]);
}

////////////////////////////////////////////////////////////////////////////////
// LodePNG decode-encode: decodes the image, then encodes it again
// colorType, bitDepth and interlace will be the same, no matter what it was
// Basically, all known chunks are interpreted and then regenerated, while all
// unknown chunks are simply placed back literally
////////////////////////////////////////////////////////////////////////////////

#include "lodepng.h"

void example8(int argc, char *argv[])
{
  std::vector<unsigned char> image;
  std::vector<unsigned char> buffer;
  LodePNG::Encoder encoder;
  LodePNG::Decoder decoder;
  
  //check if user gave a filename
  if(argc < 3)
  {
    std::cout << "please provide out and in filename" << std::endl;
    return;
  }
  
  decoder.getSettings().color_convert = 0;
  decoder.getSettings().rememberUnknownChunks = 1; //make it reproduce even unknown chunks in the saved image
  
  LodePNG::loadFile(buffer, argv[2]);
  decoder.decode(image, buffer);
  if(decoder.hasError())
  {
    std::cout << "decoder error " << decoder.getError() << ": " << LodePNG_error_text(decoder.getError()) << std::endl;
    return;
  }
  
  buffer.clear();
  
  encoder.swapInfoPng(decoder.getInfoPng());
  encoder.setInfoRaw(decoder.getInfoRaw()); //the decoder has written the PNG colortype in the infoRaw too
  encoder.settings.text_compression = 1;
  
  LodePNG_create_chunk(&encoder.getInfoPng().unknown_chunks.data[0], &encoder.getInfoPng().unknown_chunks.datasize[0], 0, "aaAa", 0);
  LodePNG_create_chunk(&encoder.getInfoPng().unknown_chunks.data[1], &encoder.getInfoPng().unknown_chunks.datasize[1], 0, "bbBb", 0);
  LodePNG_create_chunk(&encoder.getInfoPng().unknown_chunks.data[2], &encoder.getInfoPng().unknown_chunks.datasize[2], 0, "ooOo", 0);
  
  encoder.encode(buffer, image, encoder.getInfoPng().width, encoder.getInfoPng().height); /*encoder width and height used since we swapped it with that of the decoder*/
  if(encoder.hasError())
  {
    std::cout << "encoder error " << encoder.getError() << ": " << LodePNG_error_text(encoder.getError()) << std::endl;
    return;
  }
  
  LodePNG::saveFile(buffer, argv[1]);
}


////////////////////////////////////////////////////////////////////////////////
// LodePNG chunk lister: lists the types of all chunks in the PNG image
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "lodepng.h"

void example9(int argc, char *argv[]) /*list the chunks*/
{
  const char* filename = argc > 1 ? argv[1] : "test.png";
  unsigned char *chunk, *begin, *end;
  std::vector<unsigned char> buffer;
  
  LodePNG::loadFile(buffer, filename); //load the image file with given filename
  end = &buffer.back() + 1;
  begin = chunk = &buffer.front() + 8;
  
  std::cout << "type / length" << std::endl;
  std::cout << "-------------" << std::endl;
  while(chunk + 8 < end && chunk >= begin)
  {
    char type[5];
    LodePNG_chunk_type(type, chunk);
    std::cout << type << " / " << LodePNG_chunk_length(chunk) << std::endl;
    chunk = LodePNG_chunk_next(chunk);
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  //NOTE: many of these examples overwrite a file without asking so watch out what you give as parameter to the program!
  
  example1(argc, argv); //console example
  //example2(argc, argv); //SDL example (enable the #define above)
  //example3(argc, argv); //encoder example
  //example4(argc, argv); //greyscale example
  //example5(argc, argv); //palette example
  //example6(argc, argv); //ghostify example
  //example7(argc, argv); //4-bit palette example
  //example8(argc, argv); //re-encode example
  //example9(argc, argv); //chunk listing example
}
