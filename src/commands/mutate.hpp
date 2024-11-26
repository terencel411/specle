#ifndef SPECKLE_COMMANDS_MUTATE_HPP
#define SPECKLE_COMMANDS_MUTATE_HPP

#include <list>
#include <string>

class TokenStream;

void mutateParabolaCommand(TokenStream& stack);
void mutateGridNearestNeighbourCommand(TokenStream& stack);
void mutatePropagateCommand(TokenStream& stack);
void mutateFocalFarField(TokenStream& stack);
void mutateOmToTime(TokenStream& stack);

#endif // SPECKLE_COMMANDS_MUTATE_HPP