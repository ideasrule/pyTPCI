/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "module.h"


void module_list::zero() const
{
	for (vector<module *>::const_iterator it = m_l.begin(); it != m_l.end(); ++it)
	{
		(*it)->zero();
	}
}

void module_list::comment(t_warnings& w) const
{
	for (vector<module *>::const_iterator it = m_l.begin(); it != m_l.end(); ++it)
	{
		(*it)->comment(w);
	}
}
