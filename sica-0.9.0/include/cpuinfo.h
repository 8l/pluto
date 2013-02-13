/*
 * cpuinfo.h
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 */

#ifndef CPUINFO_H_
#define CPUINFO_H_

void get_register_info(unsigned int , unsigned int , char *, char *, char *, char *, const int );

void get_cpu_name(char * , char * , char * , char * , const int );

int get_intelamd_l2cache_size(char *, const int );

int get_intelamd_l2cache_line_size(char *, const int );

int get_intelamd_cacheline_size(char *, const int );

int get_intelamd_l2cache_assiociativity(char *, const int );

void get_cpu_info(int * , int * ,int * ,int *,int * );

#endif /* CPUINFO_H_ */
