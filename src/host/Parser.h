/*
 * Parser.h
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <expat.h>
#include "Cell.h"
#include "Material.h"
#include "Lattice.h"
#include "Options.h"

#ifdef DOUBLE
#define strtok_FP_PRECISION strtok_double
#else
#define strtok_FP_PRECISION strtok_float
#endif


/**
 * Parses XML Files
 */
class Parser {
private:
	std::vector<Surface *> surfaces;
	std::vector<Cell *> cells;
	std::vector<Lattice *> lattices;
	std::vector<Material *> materials;
public:
	Parser(const Options *opts);
	virtual ~Parser();
	std::vector<Surface*> getSurfaces();
	std::vector<Cell*> getCells();
	std::vector<Lattice*> getLattices();
	std::vector<Material*> getMaterials();
private:
	/* Internal functions for the XML parser, but these need to be able to
	 * access private data (and they're C functions so they can't be static
	 * class methods as they need to be passed to the C API provided by
	 * Expat) so they need to be friend'd.
	 */
	friend void XMLCALL Parser_XMLCallback_Start(void *context,
												 const XML_Char *name,
												 const XML_Char **atts);
	friend void XMLCALL Parser_XMLCallback_End(void *context,
											   const XML_Char *name);
	friend void XMLCALL Parser_XMLCallback_CData(void *context,
												 const XML_Char *s,
												 int len);
};

#endif /* POINT_H_ */
