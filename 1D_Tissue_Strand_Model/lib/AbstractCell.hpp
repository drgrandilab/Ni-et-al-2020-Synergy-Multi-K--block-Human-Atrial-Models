/*
 * AbstractCell.hpp
 * Base Class for all cell classess
 * Author: Haibo Ni
 * qiangzi.ni@gmail.com
 * Date: Fri 30 Jan 2015 14:08:11 GMT
 */

#ifndef ABSTRACTCELL_HPP
#define ABSTRACTCELL_HPP
#include "EnumSimulationCtrl.hpp"

class AbstractCell
{
public:
    AbstractCell() {};
    ~AbstractCell() {};
    TypeCell mcell_type;
    virtual void InitialiseStates(TypeCell usr_celltype);
    virtual void SetPotential(double value);
    virtual void RunSingleTimeStep(double dt, double Istim);
    virtual void ReturnPotential(double target);
};

#endif