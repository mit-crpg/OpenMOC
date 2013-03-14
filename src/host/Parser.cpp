/*
 * Parser.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 */

#include "Parser.h"
#include <vector>
#include <cstring>
#include <string>
#include <stdexcept>
#include <assert.h>

/* Verbose debugging of the parser, but doesn't follow log* format */
//#define DEBUG

/* These should really be static, but thanks to C++ that's impossible.  I've
 * still declared them up here though, as at least they can be made private!
 */
void XMLCALL Parser_XMLCallback_Start(void *context,
									  const XML_Char *name,
									  const XML_Char **atts);
void XMLCALL Parser_XMLCallback_End(void *context,
									const XML_Char *name);
void XMLCALL Parser_XMLCallback_CData(void *context,
									  const XML_Char *s,
									  int len);
/**
 * Frame to keep track of all nodes in parser
 */
enum frame_type {
	NODE_TYPE_NONE,
	NODE_TYPE_GEOMETRY,
	NODE_TYPE_MIN = NODE_TYPE_GEOMETRY,
	NODE_TYPE_CELL,
	NODE_TYPE_LATTICE,
	NODE_TYPE_TYPE,
	NODE_TYPE_DIMENSION,
	NODE_TYPE_WIDTH,
	NODE_TYPE_UNIVERSES,
	NODE_TYPE_SURFACE,
	NODE_TYPE_MATERIALS,
	NODE_TYPE_MATERIAL,
	NODE_TYPE_MAX = NODE_TYPE_MATERIAL /* Keep this in sync */
};

/**
 * Frame to keep track of nodes specific to geometry
 */
struct frame_geometry {
};

/**
 * Frame to keep track of nodes specific to materials
 */
struct frame_materials {
};

/**
 * Frame to keep track of data specific to cell
 */
struct frame_cell {
	bool has_id;
	short int id;

	bool has_fill;
	short int fill;

	bool has_material;
	short int material;

	bool has_universe;
	short int universe;

	bool has_rings;
	short int rings;

	bool has_sectors;
	short int sectors;

	int surfaces_count;
	short int *surfaces;
};

/**
 * Frame to keep track of lattice type
 */
struct frame_ttype {
	char *data;
};

/**
 * Frame to keep track of lattice dimension
 */
struct frame_dimension {
	char *data;
};

/**
 * Frame to keep track of lattice origin
 */
struct frame_origin {
	char *data;
};

/**
 * Frame to keep track of lattice width
 */
struct frame_width {
	char *data;
};

/**
 * Frame to keep track of lattice universes
 */
struct frame_universes {
	char *data;
};

/**
 * Frame to keep track of lattice information
 */
struct frame_lattice {
	bool has_id;
	short int id;

	char *type;

	short int *dimmensions;
	int dimmensions_count;

	short int *universes;
	int universes_count;

	double *width;
	int width_count;
};

/**
 * Frame to keep track of surface information
 */
struct frame_surface {
	bool has_id;
	short int id;

	char *type;

	double *coeffs;
	int coeffs_count;

	boundaryType boundary;
};

/**
 * Frame to keep track of material information
 */
struct frame_material {
	bool has_id;
	short int id;

	FP_PRECISION *sigma_a;
	int sigma_a_cnt;

	FP_PRECISION *sigma_t;
	int sigma_t_cnt;
	
	FP_PRECISION *sigma_f;
	int sigma_f_cnt;

	FP_PRECISION *nu_sigma_f;
	int nu_sigma_f_cnt;

	FP_PRECISION *chi;
	int chi_cnt;
	
	FP_PRECISION *sigma_s;
	int sigma_s_cnt;
};

/**
 * High-level frame combining all other frames
 */
struct frame {
	struct frame *parent;
	enum frame_type type;
	unsigned int depth;

	union {
		struct frame_geometry geometry;
		struct frame_materials materials;
		struct frame_cell cell;
		struct frame_lattice lattice;
		struct frame_ttype ttype;
		struct frame_dimension dimmension;
		struct frame_width width;
		struct frame_universes universes;
		struct frame_surface surface;
		struct frame_material material;
	};
};

/**
 * Stack for book-keeping inside of parser
 */
struct stack {
	struct frame *top;
	Parser *parser;
};

static inline const char *frame_type_string(enum frame_type type);

static inline struct frame *stack_push(struct stack *s, enum frame_type type);
static inline struct frame *stack_pop(struct stack *s);
static inline void stack_print(struct frame *f);
static inline void stack_print_help(struct frame *f);

static inline short int *strtok_int(const char *str, int *count);
static inline double *strtok_double(const char *str, int *count);
static inline float *strtok_float(const char *str, int *count);

static inline char *astrncat(char *orig, char *next, int len);

/**
 * Default Parser constructor
 * @param options a pointer to the options, e.g., pathes to input files
 */
Parser::Parser (const Options *opts) {
	FILE* geofile;
	FILE* matfile;
	XML_Parser parser;
	struct stack stack;
	char c;

	/* Sets up the parser */
	stack.top = NULL;
	stack.parser = this;
	parser = XML_ParserCreate(NULL); /* NULL -> system encoding */
	XML_SetUserData(parser, &stack);
	XML_SetStartElementHandler(parser, &Parser_XMLCallback_Start);
	XML_SetEndElementHandler(parser, &Parser_XMLCallback_End);
	XML_SetCharacterDataHandler(parser, &Parser_XMLCallback_CData);

	/* Assures that the input file(s) exists and is readable */
	geofile = fopen(opts->getGeometryFile(), "r");
	if (geofile == NULL) {
		log_printf(ERROR, "Given geometry file %s does not exist",
				   opts->getGeometryFile());
	}

	/* Passes single characters to the parser, which is quite slow but
	 * is the easiest for now. */
	while( EOF != (c = fgetc(geofile)) ) {
		if (XML_Parse(parser, &c, 1, false) != XML_STATUS_OK)
			log_printf(ERROR, "Expat error for geometry.xml\n");
        }

	fclose(geofile);

	/* Tells the parse we've red the end */
	XML_Parse(parser, NULL, 0, true);
	XML_ParserFree(parser);

	/* Sets up the parser */
	stack.top = NULL;
	stack.parser = this;
	parser = XML_ParserCreate(NULL); /* NULL -> system encoding */
	XML_SetUserData(parser, &stack);
	XML_SetStartElementHandler(parser, &Parser_XMLCallback_Start);
	XML_SetEndElementHandler(parser, &Parser_XMLCallback_End);
	XML_SetCharacterDataHandler(parser, &Parser_XMLCallback_CData);

	/* Assures that the input file(s) exists and is readable */
	matfile = fopen(opts->getMaterialFile(), "r");
	if (matfile == NULL) {
		log_printf(ERROR, "Given material file %s does not exist",
			   opts->getMaterialFile());
	}

	/* Passes single characters to the parser, which is quite slow but
	 * is the easiest for now. */
	while( EOF != (c = fgetc(matfile)) ) {
		if (XML_Parse(parser, &c, 1, false) != XML_STATUS_OK)
			log_printf(ERROR, "Expat error for material.xml");
        }

	fclose(matfile);

	/* Tells the parse we've reached the end */
	XML_Parse(parser, NULL, 0, true);
	XML_ParserFree(parser); 
}

/**
 * Parser Destructor
 */
Parser::~Parser() {
	materials.clear();
	surfaces.clear();
	cells.clear();
	lattices.clear();
}


std::vector<Surface*> Parser::getSurfaces() {
	return surfaces;
}

std::vector<Cell*> Parser::getCells() {
	return cells;
}

std::vector<Lattice*> Parser::getLattices() {
	return lattices;
}

std::vector<Material*> Parser::getMaterials() {
	return materials;
}


/**
 * Set handler for start tags; for each tag, validate the data, then pass 
 * data to the corresponding location on the stack.
 * @param context the stack for book-keeping
 * @param name a pointer point to the attribute name
 * @param attrs a double pointer point to the attribute value
 */
void XMLCALL Parser_XMLCallback_Start(void *context,
				      const XML_Char *name,
				      const XML_Char **attrs) {
	struct stack *s;
	struct frame *f;
	enum frame_type type;
	int i;

	s = (struct stack *)context;

	/* Checks what type of node this is */
	type = NODE_TYPE_NONE;
	for (i = (int)NODE_TYPE_MIN; i <= (int)NODE_TYPE_MAX; i++)
		if (strcmp(name, frame_type_string((enum frame_type)i)) == 0)
			type = (enum frame_type)i;

	/* Ensures that we know what type the node is */
	if (type == NODE_TYPE_NONE)
		log_printf(ERROR, "Unknown node type '%s'", name);

	/* Adds our item to the stack */
	f = stack_push(s, type);

	/* Parses every attribute */
	while (*attrs != NULL) {
		char *key, *value;
		
		/* Attributes are stored as a key-value pair, there are always
		 * two of them (one for the key, one for the value) so we can
		 * safely double-increment here.
		 */
		/* FIXME: Verify that a bad input file can't cause an odd
		 *        number of attributes.
		 */
		assert(sizeof(char) == sizeof(XML_Char));
		key = (char *)*attrs;
		attrs++;
		value = (char *)*attrs;
		attrs++;
		assert(key != NULL);
		assert(value != NULL);

		/* Does some type-specific parsing for some attributes */
		switch (f->type) {
		case NODE_TYPE_NONE:
			break;
		case NODE_TYPE_GEOMETRY:
			break;
		case NODE_TYPE_MATERIALS:
			break;
		case NODE_TYPE_CELL:
			if (strcmp(key, "id") == 0) {
				if (f->cell.has_id == true)
					log_printf(ERROR, "Cell has 2 ids");

				f->cell.has_id = true;
				f->cell.id = (short int) atoi(value);
			} else if (strcmp(key, "fill") == 0) {
				if (f->cell.has_fill == true)
					log_printf(ERROR, "Cell has 2 fills");

				if (f->cell.has_material == true) {
					log_printf(ERROR,
						   "Cell has material & fill");
				}

				f->cell.has_fill = true;
				f->cell.fill = (short int) atoi(value);
				/* Set the universe to 0 if its value is greater than 1000 (the max value allowed)
				 * by default in case this is the base universe */
				if (f->cell.universe > 1000) 
				  f->cell.universe = 0;
			} else if (strcmp(key, "rings") == 0) {
				if (f->cell.has_rings == true)
					log_printf(ERROR, "Cell has 2 ring def");

				f->cell.has_rings = true;
				f->cell.rings = (short int) atoi(value);
			} else if (strcmp(key, "sectors") == 0) {
				if (f->cell.has_sectors == true)
					log_printf(ERROR, "Cell has 2 sector def");

				f->cell.has_sectors = true;
				f->cell.sectors = (short int) atoi(value);
			} else if (strcmp(key, "material") == 0) {
				if (f->cell.has_fill == true)
					log_printf(ERROR, "Has 2 material");

				if (f->cell.has_fill == true) {
					log_printf(ERROR,
						   "Has material and fill");
				}

				f->cell.has_material = true;
				f->cell.material = (short int) atoi(value);
			} else if (strcmp(key, "universe") == 0) {
				if (f->cell.has_universe == true)
					log_printf(ERROR, "Has 2 universes");

				f->cell.has_universe = true;
				f->cell.universe = (short int) atoi(value);
			} else if (strcmp(key, "surfaces") == 0) {
				if (f->cell.surfaces != NULL)
					log_printf(ERROR, "Has 2 surfaces");

				f->cell.surfaces =
					strtok_int(value,
						   &f->cell.surfaces_count);
			} else {
				log_printf(ERROR, "Unknown attribute '%s=%s'",
					   key, value);
			}
			break;
		case NODE_TYPE_LATTICE:
			if (strcmp(key, "id") == 0) {
				if (f->lattice.has_id == true)
					log_printf(ERROR, "Lattice has 2 ids");

				f->lattice.has_id = true;
				f->lattice.id = (short int) atoi(value);
			} else {
				log_printf(ERROR, "Unknown attribute '%s=%s'",
					   key, value);
			}

			break;
		case NODE_TYPE_TYPE:
			break;
		case NODE_TYPE_DIMENSION:
		break;

		case NODE_TYPE_WIDTH:
			break;
		case NODE_TYPE_UNIVERSES:
			break;
		case NODE_TYPE_MATERIAL:
			if (strcmp(key, "id") == 0) {
				if (f->material.has_id == true)
					log_printf(ERROR, "Has 2 material ids");
				
				f->material.has_id = true;
				f->material.id = (short int) atoi(value);
			} else if (strcmp(key, "sigma_t") == 0) {
				if (f->material.sigma_t != NULL)
					log_printf(ERROR, "Has 2 sigma_t");

				f->material.sigma_t =
				  strtok_FP_PRECISION(value,
							     &f->material.sigma_t_cnt);
			} else if (strcmp(key, "sigma_f") == 0) {
				if (f->material.sigma_f != NULL)
					log_printf(ERROR, "Has 2 sigma_f");

				f->material.sigma_f =
				  strtok_FP_PRECISION(value,
								   &f->material.sigma_f_cnt);
			} else if (strcmp(key, "sigma_a") == 0) {
				if (f->material.sigma_a != NULL)
					log_printf(ERROR, "Has 2 sigma_a");

				f->material.sigma_a =
				  strtok_FP_PRECISION(value,
								   &f->material.sigma_a_cnt);
			} else if (strcmp(key, "nu_sigma_f") == 0) {
				if (f->material.nu_sigma_f != NULL)
					log_printf(ERROR, "Has 2 nu_sigma_f");

				f->material.nu_sigma_f =
				  strtok_FP_PRECISION(value,
						      &f->material.nu_sigma_f_cnt);
			} else if (strcmp(key, "chi") == 0) {
				if (f->material.chi != NULL)
					log_printf(ERROR, "Has 2 chi");

				f->material.chi =
				  strtok_FP_PRECISION(value,
								   &f->material.chi_cnt);
			} else if (strcmp(key, "sigma_s") == 0) {
				if (f->material.sigma_s != NULL)
					log_printf(ERROR, "Has 2 sigma_s");

				f->material.sigma_s =
				  strtok_FP_PRECISION(value,
						      &f->material.sigma_s_cnt);
			}
			break;
		case NODE_TYPE_SURFACE:
			if (strcmp(key, "id") == 0) {
				if (f->surface.has_id == true)
					log_printf(ERROR, "Surface has 2 ids");

				f->surface.has_id = true;
				f->surface.id = (short int) atoi(value);
			} else if (strcmp(key, "type") == 0) {
				if (f->surface.type != NULL)
					log_printf(ERROR, "Has 2 types");

				f->surface.type = strdup(value);
			} else if (strcmp(key, "coeffs") == 0) {
				if (f->surface.coeffs != NULL)
					log_printf(ERROR, "Has 2 coeffs");

				f->surface.coeffs =
					strtok_double(value,
						      &f->surface.coeffs_count);
			} else if (strcmp(key, "boundary") == 0) {
				if (f->surface.boundary != BOUNDARY_NONE)
					log_printf(ERROR, "Has 2 boundaries");

				if (strcmp(value, "reflective") == 0)
					f->surface.boundary = REFLECTIVE;
				else if (strcmp(value, "vacuum") == 0)
					f->surface.boundary = VACUUM;
				else
					log_printf(ERROR, "Only supports reflective/vacuum BC's");
			} else {
				log_printf(ERROR, "Unknown attribute '%s=%s'",
					   key, value);
			}
			break;
		}
	}
}

/**
 * Set handler for end tags; process data on the stack, allocate memory for
 * new class, construct them with corresponding data from the stack, free
 * any pointer after we are done.
 * @param context the stack for book-keeping
 * @param name a pointer point to the attribute name
 */
void XMLCALL Parser_XMLCallback_End(void *context,
				    const XML_Char *name) {
	struct stack *s;
	struct frame *f, *p;
	
	s = (struct stack *)context;
	f = stack_pop(s);
	p = s->top;

	switch (f->type) {
	case NODE_TYPE_NONE:
		break;
	case NODE_TYPE_GEOMETRY:
		break;
	case NODE_TYPE_MATERIALS:
		break;
	case NODE_TYPE_CELL:
	{
		Cell *cell;

		cell = NULL;
		if (f->cell.has_fill) {
			cell = new CellFill(f->cell.id, 
								f->cell.universe,
								f->cell.surfaces_count,
								f->cell.surfaces,
								f->cell.fill);
			if (f->cell.id == 925)
				log_printf(NORMAL, "Parser added new cell with id = %d, universe id = %d, universe fill id = %d",
															f->cell.id, f->cell.universe, f->cell.fill);
		} else if (f->cell.has_material) {
			if ( (f->cell.has_rings) && !(f->cell.has_sectors)){
				cell = new CellBasic(f->cell.id,
									 f->cell.universe,
									 f->cell.surfaces_count,
									 f->cell.surfaces,
									 f->cell.material,
									 f->cell.rings,
									 0);
			} else if (f->cell.has_sectors && !(f->cell.has_rings)){
				cell = new CellBasic(f->cell.id,
									 f->cell.universe,
									 f->cell.surfaces_count,
									 f->cell.surfaces,
									 f->cell.material,
									 0,
									 f->cell.sectors);
			} else if (f->cell.has_rings && f->cell.has_sectors) {
				cell = new CellBasic(f->cell.id, 
									 f->cell.universe,
									 f->cell.surfaces_count,
									 f->cell.surfaces,
									 f->cell.material,
									 f->cell.rings,
									 f->cell.sectors);
			} else {
				cell = new CellBasic(f->cell.id, 
									 f->cell.universe,
									 f->cell.surfaces_count,
									 f->cell.surfaces,
									 f->cell.material,
									 0, 0);		     
			} 
		} else {
			log_printf(ERROR, "Cell without material or fill");
		}
		
		if (cell != NULL)
			s->parser->cells.push_back(cell);
		else
			log_printf(ERROR, "Unknown cell type");
		
		if (f->cell.surfaces != NULL)
			free(f->cell.surfaces);
		break;
	}
	case NODE_TYPE_LATTICE:
		Lattice *lattice;

		lattice = NULL;
		if (f->lattice.has_id != true)
			log_printf(ERROR, "Lattice without id");
		if (f->lattice.dimmensions_count != 2)
			log_printf(ERROR, "Lattice without exactly 2 dimms");
		if (f->lattice.width_count != 2)
			log_printf(ERROR, "Lattice without exactly 2 widths");
		if (f->lattice.universes == NULL)
			log_printf(ERROR, "Lattice without universes");

		lattice = new Lattice((short int) f->lattice.id,
				      (short int) f->lattice.dimmensions[0],
				      (short int) f->lattice.dimmensions[1],
				      f->lattice.width[0],
				      f->lattice.width[1],
				      (short int) f->lattice.universes_count,
				      (short int*) f->lattice.universes);

		s->parser->lattices.push_back(lattice);

		if (f->lattice.type != NULL)
			free(f->lattice.type);
		if (f->lattice.dimmensions != NULL)
			free(f->lattice.dimmensions);
		if (f->lattice.universes != NULL)
			free(f->lattice.universes);
		if (f->lattice.width != NULL)
			free(f->lattice.width);
		break;
	case NODE_TYPE_TYPE:
		switch (p->type) {
		case NODE_TYPE_LATTICE:
			p->lattice.type = f->ttype.data;
			break;
		case NODE_TYPE_NONE:
		case NODE_TYPE_GEOMETRY:
		case NODE_TYPE_CELL:
		case NODE_TYPE_TYPE:	
		case NODE_TYPE_DIMENSION:
		case NODE_TYPE_WIDTH:
		case NODE_TYPE_UNIVERSES:
		case NODE_TYPE_SURFACE:
		case NODE_TYPE_MATERIAL:
		case NODE_TYPE_MATERIALS:
			log_printf(ERROR, "Unexpected type subfield");
		}
		break;
	case NODE_TYPE_DIMENSION:
		switch (p->type) {
		case NODE_TYPE_LATTICE:
			if (p->lattice.dimmensions != NULL)
				log_printf(ERROR, "Has 2 dimmensions");

			p->lattice.dimmensions =
				strtok_int(f->dimmension.data,
					   &p->lattice.dimmensions_count);
			
			free(f->dimmension.data);
			break;
		case NODE_TYPE_NONE:
		case NODE_TYPE_GEOMETRY:
		case NODE_TYPE_CELL:
		case NODE_TYPE_TYPE:	
		case NODE_TYPE_DIMENSION:
		case NODE_TYPE_WIDTH:
		case NODE_TYPE_UNIVERSES:
		case NODE_TYPE_SURFACE:
		case NODE_TYPE_MATERIALS:
		case NODE_TYPE_MATERIAL:
			log_printf(ERROR, "Unexpected dimmension subfield");
		}
		break;

		case NODE_TYPE_WIDTH:
		switch (p->type) {
		case NODE_TYPE_LATTICE:
			if (p->lattice.width != NULL)
				log_printf(ERROR, "Has 2 widths");

			p->lattice.width =
				strtok_double(f->width.data,
					      &p->lattice.width_count);
			
			free(f->width.data);
			break;
		case NODE_TYPE_NONE:
		case NODE_TYPE_GEOMETRY:
		case NODE_TYPE_CELL:
		case NODE_TYPE_TYPE:	
		case NODE_TYPE_DIMENSION:
		case NODE_TYPE_WIDTH:
		case NODE_TYPE_UNIVERSES:
		case NODE_TYPE_SURFACE:
		case NODE_TYPE_MATERIAL:
		case NODE_TYPE_MATERIALS:
			log_printf(ERROR, "Unexpected dimmension subfield");
		}
		break;
	case NODE_TYPE_UNIVERSES:
		switch (p->type) {
		case NODE_TYPE_LATTICE:
			if (p->lattice.universes != NULL)
				log_printf(ERROR, "Has 2 universes");

			p->lattice.universes =
				strtok_int(f->universes.data,
					   &p->lattice.universes_count);
			
			free(f->universes.data);
			break;
		case NODE_TYPE_NONE:
		case NODE_TYPE_GEOMETRY:
		case NODE_TYPE_CELL:
		case NODE_TYPE_TYPE:	
		case NODE_TYPE_DIMENSION:
		case NODE_TYPE_WIDTH:
		case NODE_TYPE_UNIVERSES:
		case NODE_TYPE_SURFACE:
		case NODE_TYPE_MATERIAL:
		case NODE_TYPE_MATERIALS:
			log_printf(ERROR, "Unexpected universes subfield");
		}
		break;
	case NODE_TYPE_SURFACE:
	{
		Surface *surface;

		surface = NULL;
		if (strcmp(f->surface.type, "plane") == 0) {
			if (f->surface.coeffs_count != 4)
				log_printf(ERROR, "Wrong number of coeffs");

			surface = new Plane(f->surface.id, 
					    f->surface.boundary,
					    f->surface.coeffs[0],
					    f->surface.coeffs[1],
					    f->surface.coeffs[2],
					    f->surface.coeffs[3]);
		} else if (strcmp(f->surface.type, "x-plane") == 0) {
			if (f->surface.coeffs_count != 1)
				log_printf(ERROR, "Wrong number of coeffs");

			surface = new XPlane(f->surface.id, 		       
					     f->surface.boundary,
					     f->surface.coeffs[0]);
		} else if (strcmp(f->surface.type, "y-plane") == 0) {
			if (f->surface.coeffs_count != 1)
				log_printf(ERROR, "Wrong number of coeffs");

			surface = new YPlane(f->surface.id, 
					     f->surface.boundary,
					     f->surface.coeffs[0]);
		} else if (strcmp(f->surface.type, "z-plane") == 0) {
			if (f->surface.coeffs_count != 1)
				log_printf(ERROR, "Wrong number of coeffs");

			surface = new ZPlane(f->surface.id,
					     f->surface.boundary,
					     f->surface.coeffs[0]);
		} else if (strcmp(f->surface.type, "circle") == 0) {
			if (f->surface.coeffs_count != 3)
				log_printf(ERROR, "Wrong number of coeffs");

			surface = new Circle(f->surface.id, 
					     f->surface.boundary,
					     f->surface.coeffs[0],
					     f->surface.coeffs[1],
					     f->surface.coeffs[2]);
		}
		
		if (surface != NULL)
			s->parser->surfaces.push_back(surface);
		else {
			log_printf(ERROR, "Unknown surface type '%s'",
				   f->surface.type);
		}

		if (f->surface.type != NULL)
			free(f->surface.type);
		if (f->surface.coeffs != NULL)
			free(f->surface.coeffs);
		break;
	}
	case NODE_TYPE_MATERIAL:
	{
		Material *material;

		material = new Material(f->material.id,
								f->material.sigma_a,
								f->material.sigma_a_cnt,
								f->material.sigma_t,
								f->material.sigma_t_cnt,
								f->material.sigma_f,
								f->material.sigma_f_cnt,
								f->material.nu_sigma_f,
								f->material.nu_sigma_f_cnt,
								f->material.chi,
								f->material.chi_cnt,
								f->material.sigma_s,
								f->material.sigma_s_cnt);
		
		if (material != NULL)
			s->parser->materials.push_back(material);
		else {
			log_printf(ERROR, "Material unsuccessfully built");
		}
		if (f->material.sigma_a != NULL)
			free(f->material.sigma_a);
		if (f->material.sigma_t != NULL)
			free(f->material.sigma_t);
		if (f->material.sigma_f != NULL)
			free(f->material.sigma_f);
		if (f->material.nu_sigma_f != NULL)
			free(f->material.nu_sigma_f);
		if (f->material.chi != NULL)
			free(f->material.chi);
		if (f->material.sigma_s != NULL)
			free(f->material.sigma_s);
	}
	}

	free(f);
}

/**
 * Update Stack: convert Lattice's data from uncast string to their 
 * corresponding formats
 * @param context pointer to the stack for book-keeping
 * @param str_uncast pointer to the uncast string
 * @param len length of the uncast string
 */
void XMLCALL Parser_XMLCallback_CData(void *context,
				      const XML_Char *str_uncast,
				      int len) {
	struct stack *s;
	struct frame *f;
	char *str;

	str = (char *)str_uncast;
	s = (struct stack *)context;
	f = s->top;

	switch (f->type) {
	case NODE_TYPE_NONE:
		break;
	case NODE_TYPE_GEOMETRY:
		break;
	case NODE_TYPE_MATERIALS:
		break;
	case NODE_TYPE_CELL:
		break;
	case NODE_TYPE_LATTICE:
		break;
	case NODE_TYPE_TYPE:
		break;
	case NODE_TYPE_DIMENSION:
		f->dimmension.data = astrncat(f->dimmension.data, str, len);
		break;
	case NODE_TYPE_WIDTH:
		f->width.data = astrncat(f->width.data, str, len);
		break;
	case NODE_TYPE_UNIVERSES:
		f->universes.data = astrncat(f->universes.data, str, len);
		break;
	case NODE_TYPE_SURFACE:
		break;
	case NODE_TYPE_MATERIAL:
		break;
	}
}

/**
 * Return the type of a node in string
 * @param type the enum frame type
 * @return type pointer to characters that correspond to the type
 */
const char *frame_type_string(enum frame_type type) {
	switch (type) {
	case NODE_TYPE_NONE:
		return "none";
	case NODE_TYPE_GEOMETRY:
		return "geometry";
	case NODE_TYPE_MATERIALS:
		return "materials";
	case NODE_TYPE_CELL:
		return "cell";
	case NODE_TYPE_LATTICE:
		return "lattice";
	case NODE_TYPE_TYPE:
		return "type";
	case NODE_TYPE_DIMENSION:
		return "dimension";
	case NODE_TYPE_WIDTH:
		return "width";
	case NODE_TYPE_UNIVERSES:
		return "universes";
	case NODE_TYPE_SURFACE:
		return "surface";
	case NODE_TYPE_MATERIAL:
		return "material";
	}
	
	abort();
	return NULL;
}

/**
 * Push a fram onto the stack and return it: allocate a new stack frame, 
 * initialize it based on nodes, add it to the stack
 * @param stack a pointer to the stack
 * @param type an enum frame type
 * @return frame a pointer to the stack frame
 */
struct frame *stack_push(struct stack *s, enum frame_type type) {
	struct frame *f;

	/* Allocates a new stack frame (it gets added way down at the end) */
	f = (struct frame *)malloc(sizeof(*f));
	if (f == NULL)
		log_printf(ERROR, "malloc returned NULL!");

	/* Different node types get initialized differently */
	f->type = type;
	switch (type) {
	case NODE_TYPE_NONE:
		free(f);
		log_printf(ERROR, "Tried to push an unknown node type");
		break;
	case NODE_TYPE_GEOMETRY:
		break;
	case NODE_TYPE_MATERIALS:
		break;
	case NODE_TYPE_CELL:
		f->cell.has_id = false;
		f->cell.has_fill = false;
		f->cell.has_rings = false;
		f->cell.has_sectors = false;
		f->cell.has_material = false;
		f->cell.has_universe = false;
		f->cell.surfaces = NULL;
		f->cell.surfaces_count = -1;
		break;
	case NODE_TYPE_LATTICE:
		f->lattice.has_id = false;
		f->lattice.type = NULL;
		f->lattice.dimmensions = NULL;
		f->lattice.universes = NULL;
		f->lattice.width = NULL;
		break;
	case NODE_TYPE_TYPE:
		f->ttype.data = NULL;
		break;
	case NODE_TYPE_DIMENSION:
		f->dimmension.data = NULL;
		break;
	case NODE_TYPE_WIDTH:
		f->width.data = NULL;
		break;
	case NODE_TYPE_UNIVERSES:
		f->universes.data = NULL;
		break;
	case NODE_TYPE_SURFACE:
		f->surface.has_id = false;
		f->surface.type = NULL;
		f->surface.coeffs = NULL;
		f->surface.boundary = BOUNDARY_NONE;
		break;
	case NODE_TYPE_MATERIAL:
		f->material.has_id = false;
		f->material.sigma_a = NULL;
		f->material.sigma_t = NULL;
		f->material.sigma_f = NULL;
		f->material.nu_sigma_f = NULL;
		f->material.chi = NULL;
		f->material.sigma_s = NULL;
		break;
	}

	/* We always have one depth larger than our parent */
	if (s->top == NULL)
		f->depth = 0;
	else
		f->depth = s->top->depth + 1;

	/* Actually adds this to the stack */
	f->parent = s->top;
	s->top = f;
	return f;
}

/**
 * Pop the top item off the stack and return it
 * @param stack a pointer to the stack
 * @return frame a pointer to the stack frame
 */
struct frame *stack_pop(struct stack *s) {
	struct frame *f;

	f = s->top;
	if (s->top != NULL)
		s->top = s->top->parent;

	return f;
}

/**
 * Call stack_print_help(frame) to print out the context in a frame
 * @param frame pointer to a frame
 */
void stack_print(struct frame *f) {
	stack_print_help(f);
	fprintf(stderr, "\n");
}

/**
 * Iterate through a frame and print the data to standard output
 * @param frame pointer to a frame
 */
void stack_print_help(struct frame *f) {
	unsigned int i;

	if (f == NULL)
		return;

	stack_print_help(f->parent);
	
	for (i = 0; i < f->depth; i++)
		fprintf(stderr, " ");
	fprintf(stderr, "%s", frame_type_string(f->type));
	
	switch (f->type) {
	case NODE_TYPE_NONE:
		break;
	case NODE_TYPE_GEOMETRY:
		break;
	case NODE_TYPE_MATERIALS:
		break;
	case NODE_TYPE_CELL:
		if (f->cell.has_id)
			fprintf(stderr, " id=\"%d\"", f->cell.id);
		if (f->cell.has_material)
			fprintf(stderr, " material=\"%d\"", f->cell.material);
		if (f->cell.has_fill)
			fprintf(stderr, " fill=\"%d\"", f->cell.fill);
		if (f->cell.has_rings)
			fprintf(stderr, " rings=\"%d\"", f->cell.rings);
		if (f->cell.has_sectors)
			fprintf(stderr, " sectors=\"%d\"", f->cell.sectors);
		if (f->cell.has_universe)
			fprintf(stderr, " universe=\"%d\"", f->cell.universe);
		if (f->cell.surfaces != NULL) {
			int i;

			fprintf(stderr, " surfaces=\"");
			for (i = 0; i < f->cell.surfaces_count; i++)
				fprintf(stderr, " %d", f->cell.surfaces[i]);
			fprintf(stderr, "\"");
		}

		break;
	case NODE_TYPE_LATTICE:
		if (f->lattice.has_id)
			fprintf(stderr, " id=\"%d\"", f->lattice.id);
		if (f->lattice.type != NULL)
			fprintf(stderr, " type=\"%s\"", f->lattice.type);
		if (f->lattice.dimmensions != NULL) {
			int i;

			fprintf(stderr, " dimmensions=\"");
			for (i = 0; i < f->lattice.dimmensions_count; i++) {
				fprintf(stderr, " %d",
					f->lattice.dimmensions[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->lattice.universes != NULL) {
			int i;

			fprintf(stderr, " universes=\"");
			for (i = 0; i < f->lattice.universes_count; i++) {
				fprintf(stderr, " %d",
					f->lattice.universes[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->lattice.width != NULL) {
			int i;

			fprintf(stderr, " width=\"");
			for (i = 0; i < f->lattice.width_count; i++) {
				fprintf(stderr, " %f",
					f->lattice.width[i]);
			}
			fprintf(stderr, "\"");
		}
		break;
	case NODE_TYPE_TYPE:
		break;
	case NODE_TYPE_DIMENSION:
		break;
	case NODE_TYPE_WIDTH:
		break;
	case NODE_TYPE_UNIVERSES:
		break;
	case NODE_TYPE_SURFACE:
		if (f->surface.has_id)
			fprintf(stderr, " id=\"%d\"", f->surface.id);
		if (f->surface.type != NULL)
			fprintf(stderr, " type=\"%s\"", f->surface.type);
		if (f->surface.coeffs != NULL) {
			int i;

			fprintf(stderr, " coeffs=\"");
			for (i = 0; i < f->surface.coeffs_count; i++) {
				fprintf(stderr, " %f",
					f->surface.coeffs[i]);
			}
			fprintf(stderr, "\"");
		}
		break;
	case NODE_TYPE_MATERIAL:
		if (f->material.has_id)
			fprintf(stderr, " id=\"%d\"", f->material.id);
		if (f->material.sigma_a != NULL) {
			int i;

			fprintf(stderr, " sigma_a=\"");
			for (i = 0; i < f->material.sigma_a_cnt; i++) {
				fprintf(stderr, " %f",
					f->material.sigma_a[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->material.sigma_t != NULL) {
			int i;

			fprintf(stderr, " sigma_t=\"");
			for (i = 0; i < f->material.sigma_t_cnt; i++) {
				fprintf(stderr, " %f",
					f->material.sigma_t[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->material.sigma_f != NULL) {
			int i;

			fprintf(stderr, " sigma_f=\"");
			for (i = 0; i < f->material.sigma_f_cnt; i++) {
				fprintf(stderr, " %f",
					f->material.sigma_f[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->material.nu_sigma_f != NULL) {
			int i;

			fprintf(stderr, " nu_sigma_f=\"");
			for (i = 0; i < f->material.nu_sigma_f_cnt; i++) {
				fprintf(stderr, " %f",
					f->material.nu_sigma_f[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->material.chi != NULL) {
			int i;

			fprintf(stderr, " chi=\"");
			for (i = 0; i < f->material.chi_cnt; i++) {
				fprintf(stderr, " %f",
					f->material.chi[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->material.sigma_s != NULL) {
			int i;

			fprintf(stderr, " sigma_s=\"");
			for (i = 0; i < f->material.sigma_s_cnt; i++) {
				fprintf(stderr, " %f",
					f->material.sigma_s[i]);
			}
			fprintf(stderr, "\"");
		}
		break;
	}

	fprintf(stderr, "\n");
}

/**
 * Convert a character string to a short interger
 * @param str pointer to character string
 * @param count pointer to the number of characters
 * @return arr pointer to the interger converted from the string
 */
short int *strtok_int(const char *str, int *count) {
	short int *arr;
	int cnt;
	char *duplicated;
	char *st_tmp;
	char *tok;
	int i;
			
	cnt = 0;
	duplicated = strdup(str);
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		cnt++;
	}
	free(duplicated);
	
	arr = (short int *) malloc(sizeof(*arr) * cnt);
	
	duplicated = strdup(str);
	i = 0;
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		assert(i < cnt);
		arr[i] = (short int)atoi(tok);
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		i++;
	}
	free(duplicated);

	*count = cnt;
	return arr;
}

/**
 * Convert a character string into double
 * @param str pointer to character string
 * @param count pointer to the number of characters
 * @return arr pointer to the double converted from the string
 */
double *strtok_double(const char *str, int *count) {
	double *arr;
	int cnt;
	char *duplicated;
	char *st_tmp;
	char *tok;
	int i;
			
	cnt = 0;
	duplicated = strdup(str);
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		cnt++;
	}
	free(duplicated);
	
	arr = (double*) malloc(sizeof(*arr) * cnt);
	
	duplicated = strdup(str);
	i = 0;
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		assert(i < cnt);
		arr[i] = (double)atof(tok);
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		i++;
	}
	free(duplicated);

	*count = cnt;
	return arr;
}


/**
 * Convert a character string into single precision float
 * @param str pointer to character string
 * @param count pointer to the number of characters
 * @return arr pointer to the float converted from the string
 */
float *strtok_float(const char *str, int *count) {
	float *arr;
	int cnt;
	char *duplicated;
	char *st_tmp;
	char *tok;
	int i;
			
	cnt = 0;
	duplicated = strdup(str);
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		cnt++;
	}
	free(duplicated);
	
	arr = (float*) malloc(sizeof(*arr) * cnt);
	
	duplicated = strdup(str);
	i = 0;
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		assert(i < cnt);
		arr[i] = (float)atof(tok);
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		i++;
	}
	free(duplicated);

	*count = cnt;
	return arr;
}


/**
 * Combine two strings
 * @param orig pointer to the original string
 * @param str pointer to the string to be added to the end of orig
 * @param len the length of the str
 * @return new_data the new string adding str to the end of orig
 */ 
char *astrncat(char *orig, char *str, int len) {
	char *new_data;

	if (orig == NULL) {
		new_data = strndup(str, len);
	} else {
		int olen;

		olen = strlen(orig);
		new_data = (char *) realloc(orig, olen + len + 1);
		new_data[olen+1] = '\0';
		strncat(new_data, str, len);
	}

	return new_data;
}
