#ifndef SPECKLE_COMMANDS_HDFIO_HPP
#define SPECKLE_COMMANDS_HDFIO_HPP

#include <list>
#include <string>

class TokenStream;

void createHdfCommand(TokenStream& stack);
void openHdfCommand(TokenStream& tokens);
void writeHdfCommand(TokenStream& stack);
void writeHdfTimelineCommand(TokenStream& stack);
void readHdfTimelineCommand(TokenStream& stack);
void closeHdfCommand(TokenStream& stack);
void closeHdfReadCommand(TokenStream& stack);

#endif // SPECKLE_COMMANDS_HDFIO_HPP