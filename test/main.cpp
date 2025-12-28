#include "GrunObject.h"
#include <print>
#include <vector>

int main ()
{
	GrunObject slab("rectangle", "Slab A", 30, 3.1, 0.15, "horizontal", "Stage 1");
	slab.addGrunItem("Excavator","8","","hour(s)","/2");
	slab.addGrunItem("Delivery - Sub Grade","1","","delivery(ies)","*0.5");
	slab.addGrunItem("Subgrade - Fines","0.05 * A","","m3","/2");
	slab.addGrunItem("Delivery - Steel","1","","delivery(ies)","*0.5");
	slab.addGrunItem("Ableflex - 10mm x 100mm, Stick Backed","2L2W","/25","roll(s)","*0.5");
	slab.addGrunItem("Dowel R12 450 HDG","2L1W@0.6","","bar(s)","/ 14");
	slab.addGrunItem("Mesh SL92", "2A","/ 12.5","mat(s)","* 0.66");
	slab.addGrunItem("Tie Wire (Blek)", "1", "", "roll(s)", "/20");
	slab.addGrunItem("Kahnkreet","V","","m3","/ 2.5");
	slab.addGrunItem("Labour - Secondary", "1", "", "hour(s)", "");
	slab.addGrunItem("Dowel R12 450 HDG","5L@0.6","","bar(s)","/ 14");
	slab.addGrunItem("N12 6000", "1.2/2L@1.2", "/6.0", "bars?", "");
	
	std::println("GrunObject's details: Name: {}\n\tLength: {}\tWidth: {}\tDepth: {}\tArea: {}",slab.getObjectName(),slab.getObjectProperty("length"),slab.getObjectProperty("width"),slab.getObjectProperty("depth"),slab.getObjectProperty("area"));
	std::println("GrunObject [{}] Item List information:\n{}", slab.getObjectName(), slab.getGrunItemListInfoAsString("%Y%m%d %H%M%S"));
	slab.calculateGrunObjectTotals();
	// std::println("GrunObject's Totals Data {}", slab.getGrunObjectTotalsInfoAsString());
}