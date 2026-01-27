#include <format>
#include <print>
#include <ranges>
#include <vector>
#include "GrunObject.h"


int main ()
{
	GrunObject slab("rectangle", "Slab A", 30, 3.1, 0.15, "horizontal", "Stage 1");

	slab.addGrunItem("Excavator","8","","","hour(s)","/2");
	slab.addGrunItem("Delivery - Sub Grade","1","","","delivery(ies)","*0.5");
	slab.addGrunItem("Subgrade - Fines","0.05 * A","","","m3","/2");
	slab.addGrunItem("Delivery - Steel","1","","","delivery(ies)","*0.5");
	slab.addGrunItem("Ableflex - 10mm x 100mm, Stick Backed","2L2W","","/25","roll(s)","*0.5");
	slab.addGrunItem("Dowel R12 450 HDG","2L1W@0.6","","","bar(s)","/ 14");
	slab.addGrunItem("Mesh SL92", "2A","","/ 12.5","mat(s)","* 0.66");
	slab.addGrunItem("Tie Wire (Blek)", "1","","","roll(s)","/20");
	slab.addGrunItem("Kahnkreet","V","","","m3","/ 2.5");
	slab.addGrunItem("Labour - Secondary", "1","","","hour(s)","");
	slab.addGrunItem("N12 6000","2L2W","","/5.4","length(s)","");
	slab.addGrunItem("Dowel R12 450 HDG","150","","","bar(s)","/ 14");

	// add additional relationship to GrunItem in index 5 (Dowel R12 450 HDG) and 10 (N12 6000)
	slab.addGrunItemRelationship(slab.getGrunItemByIndex(5), "5L@0.6", "5 gammon layers of dowel allegedly");
	slab.addGrunItemRelationship(slab.getGrunItemByIndex(10),"3.0*2L@0.4","X1: Extra Bars, Straight N12, 3.0m@0.4C");

/* this comment block demonstrates how to use the search methods to find relationships and their owning GrunItem locations in GrunObject::m_items
	// search all items, for a relationship
	std::string comSearch									= "";
	std::string relSearch									= "?L?W";
	// create a vector with all the indices from the slab's (GrunObject's) m_items property, using modern C++23 views and ranges
	std::vector<size_t>	searchIndices						= std::views::iota(0uz, slab.getTotalOfGrunItems()) | std::ranges::to<std::vector>();
	std::vector<RelationshipSearchResult>	searchResults	= slab.findRelationshipByStrings(searchIndices,relSearch,comSearch,false);

	std::print("Search results for Rel:'{}' and Com:'{}':",relSearch,comSearch);
	if (!searchResults.empty()) 
	{ 
		for (const auto result : searchResults)
			std::print(" {}", result);
	}

	// std::print("Searched m_items{}. Searched for Relationship '{}' and Comment '{}'. Results were: {}", searchIndices, relSearch, comSearch, searchResults | std::views::all);
	if (searchResults.empty()) { std::println(" no relationship/comment matches found."); }
	std::println();
*/
	
	std::println("GrunObject's details: Name: {}\n\tLength: {}\tWidth: {}\tDepth: {}\tArea: {}",slab.getObjectName(),slab.getObjectProperty("length"),slab.getObjectProperty("width"),slab.getObjectProperty("depth"),slab.getObjectProperty("area"));
	std::println("GrunObject [{}] Item List information:\n{}", slab.getObjectName(), slab.getGrunItemListInfoAsString("%Y%m%d %H%M%S"));
	slab.calculateGrunObjectTotals();
	std::println("GrunObject's Totals Data {}", slab.getGrunObjectTotalsInfoAsString());
}