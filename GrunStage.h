#pragma once
#include <optional>
#include <functional>
#include <vector>
#include "GrunObject.h"

#ifndef CLASS_NAME
#define CLASS_NAME "GrunStage"
#endif

class GrunStage
{
	private:
	
	enum class SpatialExponentValue {
		None 	= 0,
		Linear  = 1,
		Area    = 2,
		Volume  = 3
	};

	std::string					m_name		= "";							// the name of this GrunStage
	std::vector<GrunObject>		m_objects;									// a vector of GrunObjects this GrunStage pwns
	SpatialExponentValue		m_rateUnit	= SpatialExponentValue::None;	// the spatial unit that the Stage's cost rate uses (the SpatialExponentValue will need to be converted into a locale-specific unit of measure like m2 vs sq feet etc. depending on end-user's locale.)
	size_t						m_rateCost	= 0;							// the Stage's cost in locale currency (stored as size_t int - format to currency specs on output to UI)

	public:
	// ctr
	GrunStage(const std::string& stageName);

	size_t						addGrunObject(		const	GrunObject& 	gObject);
	size_t						calculateRateCost();
	size_t						createGrunObject(	const 	std::string&	shapeType,
													const 	std::string&	name,
															double			x,
															double			y,
															double			z,
													const	std::string&	areaType	= "horizontal",
													const	std::string&	stage		= "");
	std::optional
		<std::reference_wrapper
			<GrunObject>>		getGrunObject(		const 	size_t 			index);
	const
		std::string				getName();
	size_t						rmGrunObject(		const	int&			objIndex);
	void						setName(			const	std::string&	name);
	void						setRateUnit(		const	SpatialExponentValue& spatialUnit)	{ m_rateUnit = spatialUnit; }
	bool						setRateUnit(		const	size_t&			spatialUnitAsInt);
	size_t						size();
};