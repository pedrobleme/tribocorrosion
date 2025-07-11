import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import weibull_min
import time
from collections import defaultdict

def read_inp_file(filename):
    """Função otimizada para ler arquivos .inp"""
    with open(filename, 'r') as fid:
        lines = [line.strip() for line in fid if line.strip()]
    
    # Identificar seções
    node_start = next(i for i, line in enumerate(lines) if line.startswith('*Node'))
    elem_start = next(i for i, line in enumerate(lines) if line.startswith('*Element'))
    nset_start = next(i for i, line in enumerate(lines) if line.startswith('*Nset'))
    
    # Processar nós
    nodes = []
    for line in lines[node_start+1:elem_start]:
        parts = [p.strip() for p in line.split(',') if p.strip()]
        if len(parts) >= 4:
            nodes.append([float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])])
    
    # Processar elementos
    ele = []
    for line in lines[elem_start+1:nset_start]:
        parts = [p.strip() for p in line.split(',') if p.strip()]
        if len(parts) >= 9:
            ele.append([int(parts[0])] + [int(p) for p in parts[1:9]])
    
    return np.array(nodes), np.array(ele)

def calculate_centroids(ele, nodes):
    """Calcula centróides com fallback para lógica do MATLAB quando necessário"""
    node_dict = {int(node[0]): node[1:4] for node in nodes}
    centroids = np.zeros((len(ele), 3))
    
    for i, elem in enumerate(ele):
        node_coords = np.array([node_dict[nid] for nid in elem[1:9]])
        
        # Primeiro tenta o cálculo por média
        centroid = np.mean(node_coords, axis=0)
        
        # Verificação se o centróide está dentro dos limites do elemento
        if not (np.all(centroid >= np.min(node_coords, axis=0)) and \
               np.all(centroid <= np.max(node_coords, axis=0))):
            
            # Fallback: usa a lógica condicional do MATLAB
            dx = node_coords[0,0] - node_coords[:,0]
            dy = node_coords[0,1] - node_coords[:,1]
            dz = node_coords[0,2] - node_coords[:,2]
            
            centroid = [
                node_coords[0,0] + (np.max(dx) if np.max(dx) > 0 else np.min(dx))/2,
                node_coords[0,1] + (np.max(dy) if np.max(dy) > 0 else np.min(dy))/2,
                node_coords[0,2] + (np.max(dz) if np.max(dz) > 0 else np.min(dz))/2
            ]
        
        centroids[i] = centroid
    
    return centroids[:, 0], centroids[:, 1], centroids[:, 2]
def find_neighbors(ele):
    """Encontra elementos vizinhos compartilhando pelo menos 4 nós"""
    node_to_elems = defaultdict(list)
    for i, elem in enumerate(ele):
        for node in elem[1:9]:
            node_to_elems[node].append(i)
    
    neighbors = [[] for _ in range(len(ele))]
    for i, elem in enumerate(ele):
        neighbor_counts = defaultdict(int)
        for node in elem[1:9]:
            for neighbor_elem in node_to_elems[node]:
                if neighbor_elem != i:
                    neighbor_counts[neighbor_elem] += 1
        
        valid_neighbors = [ele[k,0] for k, count in neighbor_counts.items() if count >= 4]
        neighbors[i] = sorted(valid_neighbors)[:6]  # Limita a 6 vizinhos
    
    return neighbors

def read_volume_file(filename):
    """Lê o arquivo de volumes"""
    volumes = []
    with open(filename, 'r') as fid:
        start_reading = False
        for line in fid:
            if 'Element' in line and 'EVOL' in line:
                start_reading = True
                next(fid)
                next(fid)
                continue
            if start_reading:
                parts = line.strip().split()
                if len(parts) == 2 and parts[0].isdigit():
                    volumes.append(float(parts[1]))
    return np.array(volumes)

# Execução principal
def main():
    start_time = time.time()
    
    # 1. Carregar dados
    print("Carregando dados...")
    nodes, ele = read_inp_file('cdp.inp')
    print(f"Nós carregados: {nodes.shape}")
    print(f"Elementos carregados: {ele.shape}")
    
    # 2. Calcular centróides
    print("Calculando centróides...")
    coordx, coordy, coordz = calculate_centroids(ele, nodes)
    
    # 3. Encontrar vizinhos
    print("Encontrando vizinhos...")
    neighbors = find_neighbors(ele)
    
    # 4. Criar matriz de conectividade
    max_neighbors = 6
    con_list = np.zeros((len(ele), max_neighbors + 1))
    con_list[:, 0] = ele[:, 0]
    for i, neigh in enumerate(neighbors):
        con_list[i, 1:len(neigh)+1] = neigh
    
    # 5. Identificar superfície
    vec = np.array([1 if len(neigh) < 6 else 0 for neigh in neighbors])
    con_list_final = np.column_stack((con_list, vec))
    
    # 6. Ler volumes
    print("Lendo volumes...")
    volume_values = read_volume_file('volume.txt')
    
    # 7. Gerar corrosão
    print("Gerando valores de corrosão...")
    shape_param = 4.0
    scale_param = 10.0
    pitting_values = weibull_min.rvs(shape_param, scale=scale_param, size=len(ele))
    print(f"Mean Pitting Value: {np.mean(pitting_values):.2f}")
    print(f"Standard Deviation: {np.std(pitting_values):.2f}")
    
    # 8. Preparar dados finais
    max_cor_param = np.max(pitting_values)
    col9 = pitting_values * con_list_final[:, 7] / max_cor_param
    
    final_data = np.column_stack((
        con_list_final[:, :8].astype(int),  # 8 primeiras colunas como inteiros
        col9,                               # Valores normalizados
        volume_values,                      # Volumes
        coordx, coordy, coordz,             # Coordenadas
        pitting_values                      # Valores de corrosão originais
    ))
    
    # 9. Salvar resultados
    print("Salvando resultados...")
    np.savetxt('SimplesJob6.txt', final_data, 
               delimiter=',',
               fmt=['%d']*8 + ['%.5f']*6)
    
    print(f"Processo concluído em {time.time() - start_time:.2f} segundos")

if __name__ == "__main__":
    main()