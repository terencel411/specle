#ifndef CAV_PARSER_HPP
#define CAV_PARSER_HPP

#include <map>
#include <list>
#include <string>
#include <functional>
#include <fstream>
#include <sstream>

class TokenStream {
  protected:
    std::list<std::string> tokens;
    std::list<std::string>::iterator it;
  public:
    TokenStream(std::list<std::string> &tokens);
    virtual ~TokenStream() {}

    std::string next();
    ptrdiff_t nextInt();
    double nextDouble();
    bool hasNext();
};

/**
 * A simple parser that reads a file and executes commands based on the contents of the file.
 */
class Parser {
    public:
        typedef std::map<std::string, std::function<void(TokenStream&)>> CommandMap;
        typedef std::map<std::string, std::function<bool(TokenStream&)>> LoopCommandMap;
    private:
        CommandMap commands;
        LoopCommandMap loopCommands;

        class TockenStreamWithLoops : public TokenStream {
            private:
                std::list<std::list<std::string>::iterator> blockLabels;
            public:
                TockenStreamWithLoops(std::list<std::string> &tokens);
                void pushBlockLabel();
                void popBlockLabel();
                void loop();
        };
    public:
        Parser();
        
        /**
         * Add a command to the parser
         */
        void addCommand(std::string name, std::function<void(TokenStream&)> func);
        
        /**
         * Add a command to the parser
         */
        void addLoopCommand(std::string name, std::function<bool(TokenStream&)> func);

        /**
         * Parse the instructions from a file
         * 
         * The file contains a series of instructions in reverse polish notation.
         * Tokens are seperated by whitespace.
         * Constants are prepended with a "'" character. This means string constants may not contain whitespace.
         * 
         */
        void parse(std::string filename);
};

#endif