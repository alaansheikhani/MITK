/*=========================================================================
 
Program:   openCherry Platform
Language:  C++
Date:      $Date$
Version:   $Revision$
 
Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.
 
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.
 
=========================================================================*/

#ifndef __CHERRY_REFERENCE_EXPRESSION_H__
#define __CHERRY_REFERENCE_EXPRESSION_H__

#include "cherryDefinitionRegistry.h"

#include "Poco/DOM/Element.h"

namespace cherry {

/**
 * This class makes use of the <b>org.opencherry.core.expressions.definitions</b>
 * extension point to evaluate the current context against pre-defined
 * expressions. It provides core expression re-use.
 * 
 * @since 3.3
 */
class ReferenceExpression : public Expression {
  
private:
  
  // consider making this a more general extension manager
  // for now it's just part of the reference expression
  static DefinitionRegistry fgDefinitionRegistry;

	static DefinitionRegistry& GetDefinitionRegistry();

	static const std::string ATT_DEFINITION_ID;

	/**
	 * The seed for the hash code for all equals expressions.
	 */
	static const intptr_t HASH_INITIAL;

	std::string fDefinitionId;

	
public:
  
  ReferenceExpression(const std::string& definitionId);

	ReferenceExpression(IConfigurationElement* element);

	ReferenceExpression(Poco::XML::Element* element);

	EvaluationResult Evaluate(IEvaluationContext* context);

	void CollectExpressionInfo(ExpressionInfo* info);

	bool operator==(Expression& object);

	
protected:
  intptr_t ComputeHashCode();
};

} // namespace cherry

#endif // __CHERRY_REFERENCE_EXPRESSION_H__
