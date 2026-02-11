#pragma once
#include <vector>
#include <string>
#include <SQLiteCpp/SQLiteCpp.h>
#include "GrunItem.h"

#ifndef CLASS_NAME
#define CLASS_NAME "InventoryManager"
#endif

// forward delcarations (tells the compiler we use this class name and it will be included somewhere else beside here)
class GrunItem;

class InventoryManager
{
	public:
	// ctr
	InventoryManager(SQLite::Database& db) : m_db(db) {}

	GrunItem				getById(const int&			id);
	std::vector<GrunItem>	search(	const std::string&	nameQuery);				// bare basic search
	
	private:
	SQLite::Database&		m_db;

	void					hydrate(GrunItem& item, SQLite::Statement& query);
};