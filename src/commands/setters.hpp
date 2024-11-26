#ifndef SPECKLE_COMMANDS_SETTERS_HPP
#define SPECKLE_COMMANDS_SETTERS_HPP


#include <list>
#include <string>

class TokenStream;

void setSizeCommand(TokenStream& tokens);
void setOmSizeCommand(TokenStream& tokens);
void setDomainCommand(TokenStream& tokens);
void setLambdaCommand(TokenStream& tokens);
void setLambdaDomainCommand(TokenStream& tokens);

void setSGSpectrumCommand(TokenStream& tokens);
void setRWSpectrumPhaseCommand(TokenStream& tokens);

bool loopOmCommand(TokenStream& tokens);
bool loopXYCommand(TokenStream& tokens);

#endif // SPECKLE_COMMANDS_SETTERS_HPP