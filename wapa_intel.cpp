
// Thread 0 and 10 belong to the same physical core: 1 and 11, x and x + SMTDIST
#define SMTDIST 10

// --
// WAPA policy
void optimal_transport::apply(uint64_t current_interval, const tasklist_t &tasklist)
{
    // Apply only when the amount of intervals specified has passed
	if (current_interval % every != 0)
		return;
	if (current_interval < firstInterval)
    	return;
    if(current_interval == firstInterval) {

        // Num de combination of pairs 
        int size_n = tasklist.size();
        combinations = ( double (tasklist.size() * (tasklist.size() -1)) ) / 2.0;
		LOGINF("Num. of combinations: {}"_format(combinations));
		M = new double*[size_n];
		for(int pos = 0; pos < size_n; ++pos){ M[pos] = new double [size_n];}
        // ------
        for(int M_row = 0; M_row < size_n; M_row ++){
			for(int M_col = 0; M_col < size_n;  M_col++){
					M [M_row][M_col] = 0.0;
			}
		}

        // Initial configuration
        // Number of rounds to form the pairs
        int num_rondas = size_n - 1;
        // 1. Preparation: Vector of moving elements (1 to N-1)
        // 0 is the "Pivot" and is not included here.
        std::vector<int> moviles;
        for(int i = 0; i < size_n - 1; ++i) {
                moviles.push_back(i + 1);
        }
        // 2. Rounds Loop
        for (int r = 0; r < num_rondas; ++r) {
            
            // Temporary vector to store the mix for this round (flat: 0,1,2,3...)
            std::vector<int> mezcla_plana;
            mezcla_plana.reserve(size_n);
        
            // --- PIVOT PAIR ---
            // 0 always plays with the first element of the moving list
            mezcla_plana.push_back(0);
            mezcla_plana.push_back(moviles[0]);
            
            // --- REMAINING PAIRS ---
            // Match endpoints: the second with the last, the third with the second-to-last...
            for (int i = 1; i <= (size_n - 2) / 2; ++i) {
                mezcla_plana.push_back(moviles[i]);                  // Left
                mezcla_plana.push_back(moviles[moviles.size() - i]); // Right
            }
        
            // Print TO CHECK EVERYTHING IS GOING WELL
            // Prepared for mix of 8 applications
            // To print more applications: add mezcla_plana[x]
            LOGINF("Ronda creada: {} {} - {} {} - {} {} - {} {}"_format( 
                mezcla_plana[0], mezcla_plana[1], 
                mezcla_plana[2], mezcla_plana[3], 
                mezcla_plana[4], mezcla_plana[5], 
                mezcla_plana[6], mezcla_plana[7]
            ));

        
            // 3. Save into your tuple structure
            // NOTE: Adapt here depending on if N=10 or N=12. 
            // Case N=8
            mezclas.push_back(std::make_tuple(
                pid_t(mezcla_plana[0]),  pid_t(mezcla_plana[1]),
                pid_t(mezcla_plana[2]),  pid_t(mezcla_plana[3]),
                pid_t(mezcla_plana[4]),  pid_t(mezcla_plana[5]),
                pid_t(mezcla_plana[6]),  pid_t(mezcla_plana[7])
            ));
        
            /// 4. ROTATION (The key to the algorithm)
            // Move the first element to the end: [1, 2, 3] -> [2, 3, 1]
            std::rotate(moviles.begin(), moviles.begin() + 1, moviles.end());
        }

		size_comb = mezclas.size()-1; 
    } // End of collecting all pairs to gather initial data

    // Combinations
	int size_wk = tasklist.size();
	int num_comb = current_interval -1;
	// This phase is only applied at the beginning of the execution to fill the cost matrix
    // Warm-up phase
	if( size_comb >= num_comb)
	{
		const auto &primer = mezclas [num_comb];
		std::map<uint32_t, uint32_t>::iterator iterador;
		uint64_t ID; pid_t PID;
		// ---
		for (const auto &puntero : tasklist){
			const Task &task     = *puntero; ID  = task.id;   PID = task.pids[0];
			id_pid[ID] = PID;	
		}
		// Allocating the applications
        // Prepared for 8, 10, 12 applications
		// To adapt: Adding a case and allocating the new applications
        if (size_wk == 10 ){
            LOGINF("..         --__--__--        ..");
		    LOGINF("comb: [{} {} - {} {} - {} {} - {} {} - {} {}]"_format(std::get<0>(primer), std::get<1>(primer), std::get<2>(primer), std::get<3>(primer), std::get<4>(primer), std::get<5>(primer), std::get<6>(primer), std::get<7>(primer) , std::get<8>(primer), std::get<9>(primer) ));
		    LOGINF("..         --__--__--        ..");
            
		    for(iterador =id_pid.begin(); iterador != id_pid.end(); iterador++){
		    	int id_principal = iterador -> first;
		    	if(id_principal == std::get<0>(primer)){
		    		set_cpu_affinity_V2(0, iterador -> second);
                }else if (id_principal == std::get<1>(primer)){
		    		set_cpu_affinity_V2(10, iterador -> second);
		    	}else if (id_principal == std::get<2>(primer)){
		    		set_cpu_affinity_V2(1, iterador -> second);
		    	}else if (id_principal ==std::get<3>(primer)){
		    		set_cpu_affinity_V2(11, iterador -> second);
                }else if (id_principal == std::get<4>(primer)){
		    		set_cpu_affinity_V2(2, iterador -> second);
                }else if (id_principal == std::get<5>(primer)){
		    		set_cpu_affinity_V2(12, iterador -> second);
                }else if (id_principal == std::get<6>(primer)){
		    		set_cpu_affinity_V2(3, iterador -> second);
                }else if (id_principal == std::get<7>(primer)) {
		    		set_cpu_affinity_V2(13, iterador -> second);
                }else if (id_principal == std::get<8>(primer)) {
		    		set_cpu_affinity_V2(4, iterador -> second);
                }else if (id_principal == std::get<9>(primer)) {
		    		set_cpu_affinity_V2(14, iterador -> second);
                }
		    }
        }
        else if (size_wk == 12 ){
            LOGINF("..         --__--__--        ..");
		    LOGINF("comb: [{} {} - {} {} - {} {} - {} {} - {} {} - {} {}]"_format(std::get<0>(primer), std::get<1>(primer), std::get<2>(primer), std::get<3>(primer), std::get<4>(primer), std::get<5>(primer), std::get<6>(primer), std::get<7>(primer) , std::get<8>(primer), std::get<9>(primer), std::get<10>(primer), std::get<11>(primer) ));
		    LOGINF("..         --__--__--        ..");

		    for(iterador =id_pid.begin(); iterador != id_pid.end(); iterador++){
		    	int id_principal = iterador -> first;
		    	if(id_principal == std::get<0>(primer)){
		    		set_cpu_affinity_V2(0, iterador -> second);
                }else if (id_principal == std::get<1>(primer)){
		    		set_cpu_affinity_V2(10, iterador -> second);
		    	}else if (id_principal == std::get<2>(primer)){
		    		set_cpu_affinity_V2(1, iterador -> second);
		    	}else if (id_principal ==std::get<3>(primer)){
		    		set_cpu_affinity_V2(11, iterador -> second);
                }else if (id_principal == std::get<4>(primer)){
		    		set_cpu_affinity_V2(2, iterador -> second);
                }else if (id_principal == std::get<5>(primer)){
		    		set_cpu_affinity_V2(12, iterador -> second);
                }else if (id_principal == std::get<6>(primer)){
		    		set_cpu_affinity_V2(3, iterador -> second);
                }else if (id_principal == std::get<7>(primer)) {
		    		set_cpu_affinity_V2(13, iterador -> second);
                }else if (id_principal == std::get<8>(primer)) {
		    		set_cpu_affinity_V2(4, iterador -> second);
                }else if (id_principal == std::get<9>(primer)) {
		    		set_cpu_affinity_V2(14, iterador -> second);
                 }else if (id_principal == std::get<10>(primer)) {
		    		set_cpu_affinity_V2(5, iterador -> second);
                 }else if (id_principal == std::get<11>(primer)) {
		    		set_cpu_affinity_V2(15, iterador -> second);
		        }
            }


        }else if (size_wk == 8 ){
            LOGINF("..         --__--__--        ..");
		    LOGINF("comb: [{} {} {} {} {} {} {} {}]"_format(std::get<0>(primer), std::get<1>(primer), std::get<2>(primer), std::get<3>(primer), std::get<4>(primer), std::get<5>(primer), std::get<6>(primer), std::get<7>(primer)));
		    LOGINF("..         --__--__--        ..");

		    for(iterador =id_pid.begin(); iterador != id_pid.end(); iterador++){
		    	int id_principal = iterador -> first;
		    	if(id_principal == std::get<0>(primer)){
		    		set_cpu_affinity_V2(0, iterador -> second);
                }else if (id_principal == std::get<1>(primer)){
		    		set_cpu_affinity_V2(10, iterador -> second);
		    	}else if (id_principal == std::get<2>(primer)){
		    		set_cpu_affinity_V2(1, iterador -> second);
		    	}else if (id_principal ==std::get<3>(primer)){
		    		set_cpu_affinity_V2(11, iterador -> second);
                }else if (id_principal == std::get<4>(primer)){
		    		set_cpu_affinity_V2(2, iterador -> second);
                }else if (id_principal == std::get<5>(primer)){
		    		set_cpu_affinity_V2(12, iterador -> second);
                }else if (id_principal == std::get<6>(primer)){
		    		set_cpu_affinity_V2(3, iterador -> second);
                }else if (id_principal == std::get<7>(primer)) {
		    		set_cpu_affinity_V2(13, iterador -> second);
                }
		    }
        }

		// collect data for each application of the mix
    	pid_t taskPID;  
        uint64_t taskID;
		double ipcTotal = 0;
        // cpi obtaining values
		for (const auto &task_ptr : tasklist) {
			const Task &task     = *task_ptr;
			taskID = task.id; 
    	    std::string taskName = task.name; 
			taskPID = task.pids[0]; 
			//############################################################# 
            // Performance counters for the CPI
            double cycles       = task.stats[0].last("cpu_clk_unhalted.thread");
            double instructions = task.stats[0].last("inst_retired.any");
            // ----
            uint64_t cpu    = get_cpu_id(taskPID);
            
			// Saving the values
			for (const auto &Task_co_runner : tasklist) {
				const Task &co_runner = *Task_co_runner;
				uint64_t ID_co_runner = co_runner.id; 
				pid_t PID_co_runner   = co_runner.pids[0]; 
				uint64_t cpu_co_runner= get_cpu_id(PID_co_runner);

				if(cpu == cpu_co_runner + SMTDIST || cpu == cpu_co_runner - SMTDIST){
                    // Minimizing the CPI
                    double cpi = cycles / instructions; 
					M [taskID][ID_co_runner] += cpi; //CPI of the applications 
					M [ID_co_runner][taskID] += cpi; //CPI of the applications   

				}
                // To avoid trying to pair an application with itself
                // High cost for the diagonal of the matrix
                if (taskID == ID_co_runner ){ M [ID_co_runner][taskID] = 100.0; } 
            } // End for
			    double inst = 0, cycl = 0, ipc = -1;
			    assert((stats == "total") | (stats == "interval"));
			    if (stats == "total"){ 
                    inst = task.stats[0].sum("inst_retired.any"); cycl = task.stats[0].sum("cpu_clk_unhalted.thread");
			    } else if (stats == "interval"){
			    	inst = task.stats[0].last("inst_retired.any"); cycl = task.stats[0].last("cpu_clk_unhalted.thread");
			    }
			    ipc = inst / cycl;
			    ipcTotal += ipc;
    	} 
		LOGINF("comb_necesarias {} | comb_realizadas {}"_format(size_comb, num_comb));
		return;
	} // end loop
    LOGINF("Applying WAPA Policy : wapa_policy");
	LOGINF("Using {} stats"_format(stats));	

    std::vector<std::vector<double>> actualizar_M(size_wk, std::vector<double>(size_wk, 0.0));
    // --
	for (const auto &task_ptr : tasklist) {
		const Task &task        = *task_ptr;   uint64_t actual_taskID  = task.id; 
		pid_t actual_taskPID    = task.pids[0]; 

		id_pid[actual_taskID]   = actual_taskPID;
        uint64_t actual_cpu     = get_cpu_id(actual_taskPID);
        // Collect the performance counters
        double cycles       = task.stats[0].last("cpu_clk_unhalted.thread");
        double instructions = task.stats[0].last("inst_retired.any");
        // ----
        // Updating the cpi of each application/pair of applications
		for (const auto &Task_co_runner : tasklist) {
			const Task &co_runner = *Task_co_runner;
			uint64_t actual_ID_co_runner = co_runner.id; 
			pid_t actual_PID_co_runner   = co_runner.pids[0]; 
			uint64_t actual_cpu_co_runner= get_cpu_id(actual_PID_co_runner);
			// --
			if(actual_cpu == actual_cpu_co_runner + SMTDIST || actual_cpu == actual_cpu_co_runner - SMTDIST){
				LOGINF("[{} {}] -> [{} , {}]  -> {} - {}"_format(actual_taskID, actual_ID_co_runner, actual_taskPID, actual_PID_co_runner, actual_cpu, actual_cpu_co_runner));
                double actualizar_cpi = cycles / instructions;
				actualizar_M [actual_taskID][actual_ID_co_runner]  += actualizar_cpi; //CPI of the applications
				actualizar_M [actual_ID_co_runner][actual_taskID]  += actualizar_cpi; //CPI of the applications 
			}
		}
    } // End loop
    // Update matrix M 
    for (int i = 0; i < size_wk; ++i) {
        for (int j = 0; j < size_wk; ++j) {
            if (actualizar_M [i][j] != 0.0 ){ 
                M [i][j] = actualizar_M [i][j];
            }
        }
    }
    // Move matrix M into a vector<double> to use it as input
    for (int i = 0; i < size_wk; ++i) {
        for (int j = 0; j < size_wk; ++j) {
            M_vector.push_back(M[i][j]);
        }
    }
    // print M
    std::cout << "M = [";
    std::cout << std::endl;
    for (int i = 0; i < size_wk; ++i) {
        std::cout << "[";
        for (int j = 0; j < size_wk; ++j) {
            std::cout << M [i][j] << " ";
        }
        std::cout << "] ";
        std::cout << std::endl;
    }
    // ----------------
    //  Apply Optimal Transport
    std::vector<double> G = ot_square_regularization(a, b, M_vector, reg);
    // ----------------
    std::cout << "G = [";
    std::cout << std::endl;
    for (int i=0; i<size_wk; i++){
        std::cout << "[";
        for (int j=0; j<size_wk; j++){
            std::cout << std::fixed << std::setprecision(3) <<  (double) G[i*size_wk+j] << ", ";
        }
        std::cout << "] ";
        std::cout << std::endl;
    }
    std::cout << "]" << std::endl; 
    // ----
    std::set<int> set;
    std::set<int> used_pairs; 
    for (int G_row = 0; G_row < size_wk; G_row ++ ) {
        // check if G_row is in set
        if (set.find(G_row) != set.end()) {
            continue;
        }
        uint32_t pos_par = -1 ;
        for (int G_col  = 0; G_col < size_wk; G_col ++) {
            if (G_row == G_col) {
                continue;
            }
            if (set.find(G_col) != set.end()) {
               continue;
            }
            if (pos_par == -1) {
                pos_par = G_col;
            }
            // 
            pos_par = G[size_wk*G_row + G_col] > G[size_wk*G_row + pos_par] ? G_col : pos_par; 
        } 
         pairs_id[G_row] = pos_par; 
         pairs_id[pos_par] = G_row;
         set.insert(G_row);
         set.insert(pos_par);
         
    }
    // --------------------------
    // Allocation the pairs
	std::map<uint32_t, uint32_t> cpu_pid;  std::map<uint32_t, uint32_t> cpu_pid_fet;  
    std::map<uint32_t, uint32_t>::iterator iter_pairs;


	for(iter_pairs = pairs_id.begin(); iter_pairs != pairs_id.end(); iter_pairs++){
		uint32_t id_i = iter_pairs -> first;	uint32_t id_j = iter_pairs -> second;
		uint32_t pid_0 = (id_pid.find(id_i))->second;	uint32_t pid_1 = (id_pid.find(id_j))->second;
		// ---
		std::map<uint32_t, uint32_t>::iterator find_cpus_0;		std::map<uint32_t, uint32_t>::iterator find_cpus_1;
		uint32_t cpu_0 =  get_cpu_id(pid_0);	find_cpus_0 = cpu_pid.find(cpu_0); 
		uint32_t cpu_1 =  get_cpu_id(pid_1);    find_cpus_1 = cpu_pid.find(cpu_1);
		// ---
		if( (cpu_0 == cpu_1 + SMTDIST || cpu_0 == cpu_1 - SMTDIST) && (find_cpus_0 == cpu_pid.end() && find_cpus_1 == cpu_pid.end())){
			cpu_pid[cpu_0] = pid_0;  cpu_pid_fet [pid_0] = 1;   
			cpu_pid[cpu_1] = pid_1;  cpu_pid_fet [pid_1] = 1;  
			LOGINF("[{} {}] ~ [{} {}] ~ [{} {}] -> case 1: app1 and app2 in position"_format(id_i, id_j, pid_0, pid_1, cpu_0, cpu_1));
						 
		} else if (find_cpus_0 == cpu_pid.end() && (cpu_pid_fet.find(pid_0))-> second  != 1 ){
			uint32_t cpu_posible_2 = cpu_0 >= SMTDIST ? cpu_0 - SMTDIST : cpu_0 + SMTDIST;
			find_cpus_1 = cpu_pid.find(cpu_posible_2);
			if(find_cpus_1 == cpu_pid.end()){	
				cpu_pid[cpu_0] = pid_0;            cpu_pid_fet [pid_0] = 1; 	
			    cpu_pid[cpu_posible_2] = pid_1;    cpu_pid_fet [pid_1] = 1;  	
			}
			LOGINF("[{} {}] ~ [{} {}] ~ [{} {}] -> case 2: app1 in position but app2 have to move"_format(id_i, id_j, pid_0, pid_1,cpu_0, cpu_posible_2));	
		
		} else if (find_cpus_1 == cpu_pid.end() && (cpu_pid_fet.find(pid_1))-> second  != 1 ){
			uint32_t cpu_posible_3 = cpu_1 >= SMTDIST ? cpu_1 - SMTDIST : cpu_1 + SMTDIST;
			find_cpus_0 = cpu_pid.find(cpu_posible_3);
			if(find_cpus_0 == cpu_pid.end()){	
				cpu_pid[cpu_posible_3] = pid_0;  cpu_pid_fet [pid_0] = 1; 
				cpu_pid[cpu_1]         = pid_1;  cpu_pid_fet [pid_1] = 1;
			}
			LOGINF("[{} {}] ~ [{} {}] ~ [{} {}] -> case 3: app2 in position but app1 have to move"_format(id_i, id_j, pid_0, pid_1, cpu_posible_3, cpu_1 ));	

		} 
		else if ( (cpu_pid_fet.find(pid_0) )-> second  != 1 && (cpu_pid_fet.find(pid_1))-> second  != 1 ) {
            // None is allocated - both apps have to be allocated
			uint32_t cpus_it = 0;
			find_cpus_0 = cpu_pid.find(cpus_it);
			while (cpu_pid.find(cpus_it) != cpu_pid.end()){   cpus_it = cpus_it + 1 ;	}
			uint32_t cpu_new_0 = cpus_it; 

			uint32_t cpu_posible_4 = cpu_new_0 >= SMTDIST ? cpu_new_0 - SMTDIST : cpu_new_0 + SMTDIST;
			find_cpus_1 = cpu_pid.find(cpu_posible_4);
			if(find_cpus_1 == cpu_pid.end()){	
				cpu_pid [cpu_new_0]     = pid_0;   cpu_pid_fet [pid_0] = 1; 	
				cpu_pid [cpu_posible_4] = pid_1;   cpu_pid_fet [pid_1] = 1;
			}
			LOGINF("[{} {}] ~ [{} {}] ~ [{} {}] -> case 4: both apps have to move"_format(id_i, id_j, pid_0, pid_1, cpu_new_0, cpu_posible_4));	

	    }	
	} // end loop
    // Allocate applications to core
    std::map<uint32_t, uint32_t>::iterator move_pairs;
	for(move_pairs = cpu_pid.begin(); move_pairs != cpu_pid.end(); move_pairs++){
		uint32_t next_cpu = move_pairs -> first;
		uint32_t move_pid = move_pairs -> second;
        // ----
		set_cpu_affinity_V2(next_cpu, move_pid);
	}
    // ------ clear all maps
    cpu_pid.clear();
    cpu_pid_fet.clear();
    pairs_id.clear();
    M_vector.clear();
}


