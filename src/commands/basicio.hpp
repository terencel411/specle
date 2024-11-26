#ifndef SPECKLE_COMMANDS_BASICIO_HPP
#define SPECKLE_COMMANDS_BASICIO_HPP

#include <list>
#include <string>

class TokenStream;

void readAsciiCommand(TokenStream& stack);
void writeBinaryCommand(TokenStream& stack);
void writeFFTCommand(TokenStream& stack);

#endif // SPECKLE_COMMANDS_BASICIO_HPP