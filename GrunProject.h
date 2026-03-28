#pragma once
#include <vector>
#include <string>
#include "GrunGlobals.h"
#include "GrunStage.h"

#ifndef CLASS_NAME
#define CLASS_NAME "GrunProject"
#endif

class GrunProject {
public:
    std::string projectName;
    
    // The local snapshot of rates for THIS project
    LabourRates projectLabourRates;

    /**
     * @brief Constructor stamps the project with CURRENT global rates
     */
    GrunProject(std::string name) : projectName(name) {
        // Clone the Singleton's rates into this specific project
        projectLabourRates = GrunGlobals::getInstance().labourRates;
    }

    /**
     * @brief Adds a GrunStage object by reference to the m_stages member (vector)
     */
    void addStage(const GrunStage& stage)
	{
        m_stages.push_back(stage);
    }

	/**
	 * @brief Creates a new stage directly inside the project.
	 * @return A reference to the newly created stage so you can add objects to it immediately.
	 */
	GrunStage& createStage(const std::string& name)
	{
    	m_stages.emplace_back(name); 

		// Return a reference to the last element we just added
	    return m_stages.back();
	}

    const std::vector<GrunStage>& getStages() const { return m_stages; }

private:
    std::vector<GrunStage> m_stages;
};