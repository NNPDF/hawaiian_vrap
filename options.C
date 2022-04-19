/*
 * Options.C
 *
 *  Created on: 26 Mar 2010
 *      Author: daniel
 */

#include "options.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "pdf.h"   // for collider
#include "Vlumifns.h"   // for exchange
using namespace std;


singleValueOption::singleValueOption(const std::string& name,const std::string& value,const std::string& helpString): NamedOption(name,helpString), d_value(value) {};

bool singleValueOption::process(std::istream& is,std::string& message) const {
	Parser<std::string> parser;
	string answer=parser.parse(is);
	if ( answer == d_value  ){
		message = "" ;
		return true;
	} else {
		//if we get there, we didn't undersstand the answer
		message = "ERROR: The only possible value is: \'" + d_value + "\' for option " + getName() + "." ;
		return false;
	}
}

std::ostream& singleValueOption::printHelp(std::ostream& os) const {
	return os << d_helpString << std::string(" The only available value is: ") << d_value << ".";
}


bool yesOrNoOption::process(std::istream& is,std::string& message) const {
	string answer;
	is >> answer;
	if ( answer == "Yes" || answer == "On" || answer == "yes" || answer == "on" || answer == "YES" || answer == "ON"  ){
		d_valueToSet=true;
		message = "" ;
		return true;
	}
	if ( answer == "No" || answer == "Off" || answer == "no" || answer == "off" || answer == "NO" || answer == "OFF"  ){
		d_valueToSet=false;
		message = "" ;
		return true;
	}
	//if we get there, we didn't undersstand the answer
	message = "ERROR: could not understand the answer: \'" + answer + "\' for option " + getName() + "." ;
	return false;
}

yesOrNoOption::yesOrNoOption(const std::string& name,bool& valueToSet,const std::string& helpString): NamedOption(name,helpString), d_valueToSet(valueToSet) {};





void OptionsHandler::printHelp(ostream& os) const {
	map<string,option*>::const_iterator it = d_optionMap.begin();

	while (it != d_optionMap.end() ){
		os << (*it).first << ": \t";
		(*it).second->printHelp(os);
		os << "\n\n";
		++it;
	}

};

bool OptionsHandler::process(const std::string& line,std::string& message) const {
		if ( d_debug ){
			cout << "treating option line \'" << line << "\'" << endl;
		}
		stringstream ss(line);
		string optionName ;
		
		if (line.size() == 0 ){
			// this is an empty line
		} else if (line.size()> 0 && line[0] == '#'){
			// this is a comment, do nothing
			return true;
		} else {
			Parser<std::string> parser;
			optionName = parser.parse(ss);
			const map<string,option*>::const_iterator it = d_optionMap.find(optionName);

			if (it != d_optionMap.end() ){
				bool result=(*it).second->process(ss,message);
				if (result){
					d_state=success;
					return true;
				} else {
					d_state=failed;
					return false;
				}
			} else {
				d_state=unknown;
				message= "No option called " + optionName + " found.";
				return false;
			}
		}
	return true;
}



bool OptionsHandler::process_file(ifstream& file,string& message) const {
	char buffer[256];
	while (file.getline(buffer,256)){
		;
		int pos=-1;
		while ( buffer[++pos] == ' ');

		int length = strlen(buffer);

		if (pos != length){
			switch ( buffer[pos]){
			case '#' : {
			} break;
			default: {
	//			cout << "setting: " << buffer << endl;
				string setting(buffer);
				if (setting.size() != 0 ){
					bool success=process(setting,message);
					if ( !success ){
						return false;
					}
				}
			} break;
			}
		}
	}
	return true;
	
}

void OptionsHandler::add(option* opt){
	d_optionMap.insert(make_pair(opt->getName(),opt));
}


template class multipleValueOption<int>;
template class multipleValueOption<double>;
template class multipleValueOption<exchange>;
template class multipleValueOption<collider>;
template class multipleValueOption<string>;

template class ValueSettingOption<string>;
template class ValueSettingOption<double>;
template class ValueSettingOption<int>;



