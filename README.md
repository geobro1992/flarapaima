# Socio-ecological models of fisheries

### SOCIAL-ECOLOGICAL SYSTEMS
Ostrom's (2009) SES (social-ecological system) multi-scalar framework is a framework that seeks to integrate social and ecological systems across multiple scales to improve resource governance and management. The framework emphasizes the importance of understanding the interconnections between social and ecological systems, and how these interconnections shape resource governance and management.

The SES multi-scalar framework is composed of three core elements: (1) resource systems, (2) resource units, and (3) governance systems. Resource systems refer to the ecological components of a system, such as the physical and biological characteristics of a natural resource. Resource units are the social components of a system, including the individuals, households, communities, and organizations that use or manage the resource. Governance systems refer to the institutional arrangements and decision-making processes that shape how resources are governed and managed.

The SES multi-scalar framework also emphasizes the importance of understanding how these three elements interact across multiple scales. These scales can include local, regional, national, and international levels. The interactions between these three elements can shape the sustainability of resource governance and management, and can be affected by factors such as policy, politics, culture, and history.

### INTEGRATING POPULATION DYNAMICS
The SES multi-scalar framework can be coupled with population dynamics theory to better understand how social and ecological factors interact to shape the dynamics of populations and their habitats.

The resource systems component of the SES multi-scalar framework can be coupled with population dynamics theory by focusing on the ecological factors that influence population dynamics. For example, understanding the size and structure of fish populations in a fishery ecosystem, as well as their reproductive rates and interactions with other species, can help to inform sustainable management strategies for the fishery.

The resource units component of the SES multi-scalar framework can be coupled with population dynamics theory by focusing on the social factors that influence population dynamics. For example, understanding the fishing practices and behaviors of fishing communities, as well as their social and economic incentives for fishing, can help to inform policies and programs that promote sustainable fishing practices.

The governance systems component of the SES multi-scalar framework can be coupled with population dynamics theory by focusing on the institutional arrangements and decision-making processes that shape resource management. For example, understanding how fishing regulations and policies are developed, enforced, and modified over time can help to inform more effective governance strategies that promote sustainable fishing practices and healthy fish populations.


### Contents
We explore three scenarios for fisheries that experience unreported harvest, each scenario comes with scripts to 1) build the model in JAGS code, 2) run the Baysian algorithm, and 3) plot parameter estimates.

**Scenario 1** - illegal harvest is proportional to fish abundance and fishing pressure is unknown. The scripts to run and evaluate models are provided in the "_IPM/ih_proportional_" folder.\
**Scenario 2** - illegal harvest is proportional to fish abundance and fishing pressure is known. The scripts to run and evaluate models are provided in the "_IPM/multinomial_" folder.\
**Scenario 3** - illegal harvest is unrelated to fish abundance and fishing pressure in unknown. The scripts to run and evaluate models are provided in the "_IPM/ih_constant_" folder.

![Conceptual Flow Chart](chans_conceptual_fig.tif)

Figure 1: Decision tree for selecting a fisheries model. Choices are made based on the relationship between unreported harvest and fish abundance, the availability of data to precisely estimate legal fishing, and the ability of social information to estimate levels of rule compliance."
