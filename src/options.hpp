/*
 * Options.hpp
 *
 *  Created on: 26 Mar 2010
 *      Author: daniel
 */


#include <string>
#include <iostream>
#include <map>
#include <sstream>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <cassert>
#include <string.h>
#include <typeinfo>
#include <vector>


#include "options.h"


template <class T> bool CheckInput(const std::string& str);

template <> bool CheckInput<int>(const std::string& str){
	char buffer[100];
	assert(str.size()< 100 );
	strcpy(buffer,str.c_str());

	char *end_ptr = buffer+str.size();
	strtol(buffer, &end_ptr,0);
	int diff=end_ptr-&buffer[0];
	if ( diff != str.size()){
		return false;
	} else {
		return true;
	}

}

template <> bool CheckInput<double>(const std::string& str){
	char buffer[100];
	assert(str.size()< 100 );
	strcpy(buffer,str.c_str());

	char *end_ptr = buffer+str.size();
	strtod(buffer, &end_ptr);
	int diff=end_ptr-&buffer[0];
	if ( diff != str.size()){
		return false;
	} else {
		return true;
	}

}


template <class T> bool CheckInput(const std::string& str){
	return true;
}

template <class T> T Parser<T>::parse(std::istream& is) const {
	std::string str;
	if (!is){
		throw  ParserNoInput() ;
	}
	is >> str;
	if ( str.size() == 0){
		throw  ParserNoInput() ;
	}
	std::stringstream ss(str);
	T value;
	ss >> value;
	if ( ss.fail() ){
		throw  ParserError(str) ;
	} else {
		if ( !CheckInput<T>(str) ){
			throw  ParserTypeError(str) ;
		}
		return value;
	}
}


template <class TargetType> ValueSettingOption<TargetType>::ValueSettingOption(const std::string& name, TargetType& valueToSet,
		const std::string& helpString): NamedOption(name,helpString), d_valueToSet(valueToSet) {
}


template <class T> inline std::string MyClassName();

template <> inline std::string MyClassName<double>(){
	return "double";
}
template <class T> inline std::string MyClassName(){
	return typeid(T).name();
}
;


template <class TargetType> bool ValueSettingOption<TargetType>::process(std::istream& is,std::string& message) const {
	TargetType value;
	Parser<TargetType> parser;

	try {
		value = parser.parse(is);
		d_valueToSet=value ;
		message = "" ;
		return true;
	}
	catch ( ParserNoInput& pni ) {
		message = "Unexpected end of input while looking for a value for option " + getName() + ". \n" ;
		return false;
	}
	catch ( ParserError& pte ) {
		message = "Unexpected value for option " + getName() + " : \'" + pte.d_str  + "\'. Expected value of type " + MyClassName<TargetType>() +"\n" ;
		return false;
	}
}


template <class TargetType> multipleValueOption<TargetType>::multipleValueOption(const std::string& name, TargetType& valueToSet,
		const std::string& name1, const TargetType& value1,
		const std::string& helpString): NamedOption(name,helpString), d_valueToSet(valueToSet) {
	d_nameMap.emplace(name1,value1);
}
template <class TargetType> multipleValueOption<TargetType>::multipleValueOption(const std::string& name, TargetType& valueToSet,
		const std::string& name1, const TargetType& value1,
		const std::string& name2, const TargetType& value2,
		const std::string& helpString): NamedOption(name,helpString), d_valueToSet(valueToSet) {
	d_nameMap.emplace(name1,value1);
	d_nameMap.emplace(name2,value2);
}
template <class TargetType> multipleValueOption<TargetType>::multipleValueOption(const std::string& name, TargetType& valueToSet,
		const std::string& name1, const TargetType& value1,
		const std::string& name2, const TargetType& value2,
		const std::string& name3, const TargetType& value3,
		const std::string& helpString): NamedOption(name,helpString), d_valueToSet(valueToSet) {
	d_nameMap.emplace(name1,value1);
	d_nameMap.emplace(name2,value2);
	d_nameMap.emplace(name3,value3);
}
template <class TargetType> multipleValueOption<TargetType>::multipleValueOption(const std::string& name, TargetType& valueToSet,
		const std::string& name1, const TargetType& value1,
		const std::string& name2, const TargetType& value2,
		const std::string& name3, const TargetType& value3,
		const std::string& name4, const TargetType& value4,
		const std::string& helpString): NamedOption(name,helpString), d_valueToSet(valueToSet) {
	d_nameMap.emplace(name1,value1);
	d_nameMap.emplace(name2,value2);
	d_nameMap.emplace(name3,value3);
	d_nameMap.emplace(name4,value4);
}
template <class TargetType> multipleValueOption<TargetType>::multipleValueOption(const std::string& name, TargetType& valueToSet,
		const std::string& name1, const TargetType& value1,
		const std::string& name2, const TargetType& value2,
		const std::string& name3, const TargetType& value3,
		const std::string& name4, const TargetType& value4,
		const std::string& name5, const TargetType& value5,
		const std::string& helpString): NamedOption(name,helpString), d_valueToSet(valueToSet) {
	d_nameMap.emplace(name1,value1);
	d_nameMap.emplace(name2,value2);
	d_nameMap.emplace(name3,value3);
	d_nameMap.emplace(name4,value4);
	d_nameMap.emplace(name5,value5);
}
template <class TargetType> multipleValueOption<TargetType>::multipleValueOption(const std::string& name, TargetType& valueToSet,
		const std::string& name1, const TargetType& value1,
		const std::string& name2, const TargetType& value2,
		const std::string& name3, const TargetType& value3,
		const std::string& name4, const TargetType& value4,
		const std::string& name5, const TargetType& value5,
		const std::string& name6, const TargetType& value6,
		const std::string& helpString): NamedOption(name,helpString), d_valueToSet(valueToSet) {
	d_nameMap.emplace(name1,value1);
	d_nameMap.emplace(name2,value2);
	d_nameMap.emplace(name3,value3);
	d_nameMap.emplace(name4,value4);
	d_nameMap.emplace(name5,value5);
	d_nameMap.emplace(name6,value6);
}


template <class TargetType> bool multipleValueOption<TargetType>::process(std::istream& is,std::string& message) const {
	std::string name;
	Parser<std::string> parser;
	name = parser.parse(is);
	typename std::map<std::string,TargetType>::const_iterator beg = d_nameMap.begin();
	typename std::map<std::string,TargetType>::const_iterator end = d_nameMap.end();
	typename std::map<std::string,TargetType>::const_iterator it = d_nameMap.find(name);

	if  ( it != end ){
		d_valueToSet=(*it).second ;
		message = "" ;
		return true;
	}
	else {
		message = "ERROR: unsupported " + getName() + " value \'" + name + "\'. The only supported values are \n   " ;
		typename std::map<std::string,TargetType>::const_iterator it;
		for (it = beg; it != end;++it ){
			message += (*it).first ;
			message += "\n   " ;
		}
		return false;
	}
}



