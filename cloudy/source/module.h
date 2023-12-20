/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MODULE_H_
#define MODULE_H_

class module;
class t_warnings;

class module_list : public Singleton<module_list>
{
	friend class Singleton<module_list>;
public:
	vector<module *> m_l;
protected:
	module_list() {}
public:
	void use(module *m)
	{
		m_l.push_back(m);
	}
	void zero() const;
	void comment(t_warnings&) const;
};

class module
{
public:
	module()
	{
		module_list::Inst().use(this);
	}
	virtual void zero() = 0;
	virtual void comment(t_warnings&) = 0;
	virtual const char* chName() const = 0;
	virtual ~module() {}
};

#endif /* MODULE_ */
