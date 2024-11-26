#include "parse.hpp"

#include <cstddef>

TokenStream::TokenStream(std::list<std::string> &tokens) {
    this->tokens = tokens;
    it = this->tokens.begin();
}
std::string TokenStream::next() {
    return *it++;
}

ptrdiff_t TokenStream::nextInt() {
    return std::stoll(next());
}
double TokenStream::nextDouble() {
    return std::stod(next());
}

bool TokenStream::hasNext() {
    return it != tokens.end();
}

Parser::Parser() {}

void Parser::addCommand(std::string name, std::function<void(TokenStream&)> func) {
    commands[name] = func;
}

void Parser::addLoopCommand(std::string name, std::function<bool(TokenStream&)> func) {
    loopCommands[name] = func;
}

void Parser::parse(std::string filename) {
    std::ifstream infile(filename);
    std::string line;
    std::list<std::string> tokens;

    while (std::getline(infile, line)) {
        
        std::istringstream iss(line);
        std::string token;

        while (iss >> token) {
          tokens.push_back(token);
        }
    }

    TockenStreamWithLoops tokenStream(tokens);

    while (tokenStream.hasNext()) {
        std::string token = tokenStream.next();
        if (token == "do") {
            tokenStream.pushBlockLabel();
        } else if (commands.find(token) != commands.end()) {
            commands[token](tokenStream);
        } else if (loopCommands.find(token) != loopCommands.end()) {
            if (loopCommands[token](tokenStream)) {
                tokenStream.loop();
            } else {
                tokenStream.popBlockLabel();
            }
        } else {
            throw std::runtime_error("Unknown instruction: " + token);
        }
    }
}

Parser::TockenStreamWithLoops::TockenStreamWithLoops(std::list<std::string> &tokens) : TokenStream(tokens) {}

void Parser::TockenStreamWithLoops::pushBlockLabel() {
    blockLabels.push_back(it);
}

void Parser::TockenStreamWithLoops::popBlockLabel() {
    blockLabels.pop_back();
}

void Parser::TockenStreamWithLoops::loop() {
    it = blockLabels.back();
}