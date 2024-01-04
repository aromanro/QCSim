#include <iostream>
#include <iterator>

#include "Tests.h"
#include "VQE.h"


// tests values follow the Lesson 11 from
// "Fundamentals In Quantum Algorithms: A Tutorial Series Using Qiskit Continued" by Daniel Koch, Saahil Patel, Laura Wessing, Paul M. Alsing
// https://arxiv.org/abs/2008.10647

bool VQETests()
{
	std::uniform_real_distribution<> dist_zo(0., 1.);

	std::cout << "\nVQE Tests" << std::endl;

	std::cout << "\nOne qubit" << std::endl;

	VQE::VariationalQuantumEigensolver<> vqeSingleQubit(1);
	VQE::PauliStringVQE<> Xterm(1);
	Xterm.setOperatorForQubit(0, PauliString::PauliString::PauliOp::opX);
	Xterm.setCoefficient(3.0);
	VQE::PauliStringVQE<> Yterm(1);
	Yterm.setOperatorForQubit(0, PauliString::PauliString::PauliOp::opY);
	Yterm.setCoefficient(-2.0);
	VQE::PauliStringVQE<> Zterm(1);
	Zterm.setOperatorForQubit(0, PauliString::PauliString::PauliOp::opZ);
	Zterm.setCoefficient(1.0);

	vqeSingleQubit.AddTerm(Xterm);
	vqeSingleQubit.AddTerm(Yterm);
	vqeSingleQubit.AddTerm(Zterm);

	const double theta = dist_zo(gen) * M_PI / 2;
	const double phi = dist_zo(gen) * M_PI;
	const double radius = 0.35;
	
	std::vector<std::vector<double>> vertices(3);

	const double R = dist_zo(gen) * 2. * M_PI / 3.;
	for (int v = 0; v < static_cast<int>(vertices.size()); ++v)
	{
		const double angle = R + v * 2. * M_PI / 3.;
		vertices[v].push_back(theta + radius * cos(angle));
		vertices[v].push_back(phi + radius * sin(angle));
	}

	// this is close to the optimal solution, at least it should stay around that
	// result is around -3.8
	//vertices[0][0] = 1.8145;
	//vertices[0][1] = 2.5429;

	vqeSingleQubit.SetVertices(vertices);

	vqeSingleQubit.Execute();

	double lowestEnergy = vqeSingleQubit.GetMinEnergy();

	std::vector<double> vertex = vqeSingleQubit.GetMinVertex();
	std::cout << "Theta: " << vertex[0] << " Phi : " << vertex[1] << std::endl;
	std::cout << "Energy: " << lowestEnergy << " Best found should be around -3.8!" << std::endl;
	
	if (lowestEnergy > -3.)
	{
		std::cout << "Energy is too high!" << std::endl;
		return false;
	}

	std::cout << "\nTwo qubits" << std::endl;

	VQE::VariationalQuantumEigensolver<> vqeDoubleQubits(2);
	VQE::PauliStringVQE<> XYterm(2);
	XYterm.setOperatorForQubit(0, PauliString::PauliString::PauliOp::opX);
	XYterm.setOperatorForQubit(1, PauliString::PauliString::PauliOp::opY);
	XYterm.setCoefficient(3.0);
	VQE::PauliStringVQE<> ZZterm(2);
	ZZterm.setOperatorForQubit(0, PauliString::PauliString::PauliOp::opZ);
	ZZterm.setOperatorForQubit(1, PauliString::PauliString::PauliOp::opZ);
	ZZterm.setCoefficient(-2.0);

	vqeDoubleQubits.AddTerm(XYterm);
	vqeDoubleQubits.AddTerm(ZZterm);

	vqeDoubleQubits.SetTerminateLimit(15);

	vertices.clear();
	vertices.resize(9);

	std::vector<double> params;
	for (int i = 0; i < 4; ++i)
	{
		params.push_back(dist_zo(gen) * M_PI / 2);
		params.push_back(dist_zo(gen) * M_PI);
	}

	for (int v = 0; v < static_cast<int>(vertices.size()); ++v)
		for (int i = 0; i < static_cast<int>(params.size()); ++i)
		{
			const double R = (0.4 + 0.8 * dist_zo(gen)) * (dist_bool(gen) ? -1 : 1);
			vertices[v].push_back(params[i] + R);
		}
	
	vqeDoubleQubits.SetVertices(vertices);
	vqeDoubleQubits.Execute();

	lowestEnergy = vqeDoubleQubits.GetMinEnergy();

	std::cout << "Energy: " << lowestEnergy << " Best found should be around -5!" << std::endl;

	if (lowestEnergy > -4)
	{
		std::cout << "Energy is too high!" << std::endl;
		return false;
	}


	return true;
}