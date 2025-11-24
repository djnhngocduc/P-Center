# Literature Review: P-Center Problem and SAT Encoding

## 1. Giới thiệu

Bài toán P-Center là một bài toán cổ điển trong lĩnh vực Nghiên cứu Vận trù học (Operations Research) và tối ưu hóa, liên quan đến việc xác định vị trí tối ưu của các cơ sở để phục vụ một tập hợp các điểm nhu cầu. Mục tiêu chính của bài toán này là giảm thiểu khoảng cách tối đa giữa bất kỳ điểm nhu cầu nào và cơ sở gần nhất của nó. Điều này có nghĩa là bài toán P-Center tìm cách đảm bảo rằng không có điểm nhu cầu nào quá xa cơ sở phục vụ nó, từ đó tối ưu hóa sự công bằng trong phân phối dịch vụ.

### 1.1. Phát biểu bài toán

Bài toán P-Center có thể được phát biểu như sau: Cho một tập hợp các điểm nhu cầu (khách hàng) và một tập hợp các vị trí tiềm năng để đặt cơ sở, hãy xác định vị trí của P cơ sở sao cho khoảng cách tối đa giữa một điểm nhu cầu và cơ sở gần nhất của nó được giảm thiểu.

Một cách hình thức, bài toán có thể được mô hình hóa như sau:

Giả sử:
*   $N$ là tập hợp các điểm nhu cầu.
*   $M$ là tập hợp các vị trí tiềm năng để đặt cơ sở.
*   $d_{ij}$ là khoảng cách giữa điểm nhu cầu $i$ và vị trí cơ sở tiềm năng $j$.
*   $P$ là số lượng cơ sở cần đặt.

Mục tiêu là chọn một tập hợp $S \subseteq M$ gồm $P$ vị trí cơ sở sao cho:

$\min \max_{i \in N} \min_{j \in S} d_{ij}$

### 1.2. Ứng dụng

Bài toán P-Center có nhiều ứng dụng thực tế trong các lĩnh vực khác nhau, bao gồm:

*   **Dịch vụ khẩn cấp:** Xác định vị trí tối ưu cho các trạm cứu thương, trạm cứu hỏa, đồn cảnh sát để đảm bảo thời gian phản ứng nhanh nhất đến bất kỳ địa điểm nào trong khu vực phục vụ.
*   **Logistics và quản lý chuỗi cung ứng:** Lập kế hoạch vị trí kho bãi, trung tâm phân phối để giảm thiểu thời gian giao hàng tối đa đến khách hàng.
*   **Y tế:** Xác định vị trí bệnh viện, phòng khám để đảm bảo mọi bệnh nhân đều có thể tiếp cận dịch vụ y tế trong khoảng cách hợp lý.
*   **Viễn thông:** Đặt các trạm phát sóng di động (cell tower) để tối ưu hóa vùng phủ sóng và giảm thiểu khoảng cách tối đa đến người dùng.
*   **Quy hoạch đô thị:** Xác định vị trí các tiện ích công cộng như trường học, thư viện, công viên để phục vụ tốt nhất cho cộng đồng.

### 1.3. Các mô hình/formulations cơ bản

Bài toán P-Center có thể được giải quyết bằng nhiều phương pháp khác nhau, bao gồm các phương pháp chính xác và các phương pháp heuristic. Các phương pháp chính xác cung cấp giải pháp tối ưu nhưng có thể tốn kém về mặt tính toán đối với các trường hợp lớn. Các phương pháp heuristic, mặt khác, cung cấp các giải pháp gần tối ưu trong thời gian tính toán hợp lý.

#### 1.3.1. Mô hình lập trình nguyên (Integer Programming - IP)

Bài toán P-Center có thể được xây dựng thành một mô hình Lập trình nguyên (IP) như sau:

**Biến quyết định:**
*   $y_j \in \{0, 1\}$: bằng 1 nếu cơ sở $j$ được chọn, 0 nếu ngược lại.
*   $x_{ij} \in \{0, 1\}$: bằng 1 nếu điểm nhu cầu $i$ được gán cho cơ sở $j$, 0 nếu ngược lại.
*   $R$: biến liên tục biểu thị khoảng cách tối đa (bán kính).

**Mục tiêu:**

$\min R$

**Các ràng buộc:**

1.  $\sum_{j \in M} y_j = P$ (Chọn đúng $P$ cơ sở)
2.  $\sum_{j \in M} x_{ij} = 1 \quad \forall i \in N$ (Mỗi điểm nhu cầu được gán cho đúng một cơ sở)
3.  $x_{ij} \le y_j \quad \forall i \in N, j \in M$ (Chỉ có thể gán điểm nhu cầu cho cơ sở đã chọn)
4.  $d_{ij} x_{ij} \le R \quad \forall i \in N, j \in M$ (Khoảng cách từ điểm nhu cầu đến cơ sở được gán không vượt quá $R$)

#### 1.3.2. Các phương pháp giải quyết

Các phương pháp giải quyết bài toán P-Center có thể được phân loại thành:

*   **Phương pháp chính xác (Exact Methods):** Đảm bảo tìm được giải pháp tối ưu toàn cục.
    *   Lập trình nguyên (Integer Programming)
    *   Branch and Bound
    *   Cutting Plane
*   **Phương pháp heuristic và metaheuristic (Heuristic and Metaheuristic Approaches):** Tìm kiếm giải pháp tốt trong thời gian hợp lý, không đảm bảo tối ưu toàn cục.
    *   Thuật toán tham lam (Greedy Algorithm)
    *   Tìm kiếm cục bộ (Local Search)
    *   Các metaheuristic (ví dụ: Genetic Algorithm, Simulated Annealing, Tabu Search, Variable Neighborhood Search).

Trong các phần tiếp theo, chúng ta sẽ đi sâu vào các nghiên cứu cụ thể sử dụng các phương pháp này, đặc biệt tập trung vào SAT encoding và các phương pháp liên quan.

## 2. Effective Approaches to Solve P-Center Problem via Set Covering and SAT (Liu et al., 2020)

*   **Phương pháp:** Giải quyết bài toán P-Center bằng cách chuyển đổi nó thành một chuỗi các bài toán con Set Covering, sau đó mã hóa các bài toán con này sang định dạng CNF (Conjunctive Normal Form) và giải quyết bằng các bộ giải SAT (Satisfiability) hiện đại.
    *   **Quy trình:**
        1.  Chuyển đổi bài toán P-Center thành một chuỗi các bài toán con Set Covering.
        2.  Đơn giản hóa các bài toán con bằng các quy tắc giảm dữ liệu.
        3.  Mã hóa các bài toán con đã đơn giản hóa thành định dạng CNF bằng hai phương pháp mã hóa khác nhau.
        4.  Sử dụng các bộ giải SAT hiện đại để giải quyết các bài toán CNF.
*   **Input:** Bài toán P-Center trên đồ thị vô hướng.
*   **Datasets:** 70 trường hợp benchmark cổ điển của bài toán P-Center.
*   **Công cụ/Thư viện/Thuật toán:**
    *   **Bộ giải SAT:** Các bộ giải SAT hiện đại (state-of-the-art SAT solvers).
    *   **Phương pháp mã hóa:** Hai chế độ mã hóa (sequential counter encoding mode và parallel counter encoding mode) để duy trì ràng buộc cardinality.
    *   **Quy tắc giảm dữ liệu:** Các quy tắc giảm dữ liệu được đề xuất bởi Alber et al. [18].
*   **Kết quả:** Cải thiện kết quả tốt nhất đã biết trước đây cho 3 trường hợp sử dụng bộ giải SAT heuristic và chứng minh tính tối ưu cho 59 trường hợp sử dụng bộ giải SAT chính xác. Phương pháp này cho phép giải quyết các bài toán con song song.




## 3. Compact MILP formulations for the p-center problem (Ales and Elloumi, 2023)

*   **Phương pháp:** Đề xuất hai công thức MILP (Mixed-Integer Linear Programming) mới và nhỏ gọn cho bài toán p-center. Công thức đầu tiên là sự cải tiến của một công thức trước đó, giảm đáng kể số lượng ràng buộc trong khi vẫn giữ nguyên giá trị tối ưu của phần thư giãn tuyến tính. Công thức thứ hai chứa ít biến và ràng buộc hơn nhưng có ràng buộc thư giãn tuyến tính yếu hơn.
    *   **Quy trình:**
        1.  Đề xuất hai công thức MILP mới: (CP1) và (CP2).
        2.  Giới thiệu một thuật toán hai bước để giải quyết bài toán p-center hiệu quả hơn, bao gồm việc tính toán các cận dưới và cận trên mạnh và giảm đáng kể kích thước của các công thức.
*   **Input:** Bài toán p-center với N khách hàng và M vị trí cơ sở tiềm năng.
*   **Datasets:** Các trường hợp từ OR-Library.
*   **Công cụ/Thư viện/Thuật toán:**
    *   **Bộ giải:** CPLEX 12.7 (Java API).
    *   **Phần cứng:** Intel XEON E3-1280 với 3.5 GHz và 32GB RAM.
    *   **Thuật toán:** Thuật toán hai bước để giải quyết bài toán p-center.
*   **Kết quả:** Công thức (CP1) cho thấy hiệu suất tốt nhất về mọi mặt (số lượng biến, ràng buộc, thời gian giải quyết) so với các công thức khác, bao gồm cả (P1) và (P2). Công thức (CP2) nhỏ gọn nhất nhưng có ràng buộc LP yếu hơn.




## 4. A CP/LS Heuristic Method for Maxmin and Minmax Location Problems with Distance Constraints (Iosif et al., 2024)

*   **Phương pháp:** Nghiên cứu một biến thể của bài toán p-dispersion và p-center có ràng buộc khoảng cách giữa các cơ sở. Đề xuất một bộ giải CP (Constraint Programming) không hoàn chỉnh sử dụng heuristic tham lam để cắt tỉa các nhánh, kết hợp với tìm kiếm cục bộ (Local Search) để ước tính cận tốt hơn tại mỗi nút.
    *   **Quy trình:**
        1.  Sử dụng heuristic tham lam để cắt tỉa các nhánh trong bộ giải CP.
        2.  Áp dụng tìm kiếm cục bộ (LS) để có được ước tính tốt hơn về cận tại mỗi nút, dẫn đến việc cắt tỉa tập trung hơn.
        3.  So sánh khung phương pháp này với các mô hình ILP và CP được triển khai trong Gurobi và OR-Tools.
*   **Input:** Bài toán p-center với ràng buộc khoảng cách, bao gồm các vị trí khách hàng và cơ sở tiềm năng được đặt ngẫu nhiên trên lưới, hoặc dựa trên các trường hợp benchmark p-median.
*   **Datasets:** Các trường hợp được tạo ngẫu nhiên trên lưới và các trường hợp benchmark p-median.
*   **Công cụ/Thư viện/Thuật toán:**
    *   **Bộ giải CP:** Bộ giải CP không hoàn chỉnh tùy chỉnh.
    *   **Bộ giải ILP:** Gurobi.
    *   **Bộ giải CP:** CP-SAT OR-Tools.
    *   **Heuristic:** Heuristic tham lam, tìm kiếm cục bộ (Local Search).
*   **Kết quả:** Bộ giải CP không hoàn chỉnh cho thấy hiệu suất mạnh mẽ hơn so với Gurobi và OR-Tools trên các bài toán lớn hơn, và tìm thấy nhiều giải pháp tối ưu hoặc gần tối ưu nhanh hơn trên các bài toán dựa trên p-median.




## 5. Exploiting flat subspaces in local search for p-Center problem and two fault-tolerant variants (Mousavi, 2023)

*   **Phương pháp:** Đề xuất các thuật toán tìm kiếm cục bộ (Local Search) cho bài toán p-Center và hai biến thể chịu lỗi của nó (α-Neighbour p-Center và p-Next Center). Các thuật toán được đề xuất chia sẻ một thiết kế chung: tích hợp tìm kiếm cục bộ cải tiến đầu tiên với các chiến lược để khai thác các không gian phẳng trong không gian tìm kiếm.
*   **Input:** Bài toán p-Center và các biến thể chịu lỗi của nó.
*   **Datasets:** Sử dụng các bộ dữ liệu tiêu chuẩn.
*   **Công cụ/Thư viện/Thuật toán:**
    *   **Thuật toán:** Tìm kiếm cục bộ (Local Search), Metaheuristic.
*   **Kết quả:** Thuật toán được đề xuất cho p-Center vượt trội hơn metaheuristic hiện đại nhất cho bài toán này. Thuật toán được đề xuất cho p-Next Center cũng vượt trội hơn metaheuristic hiện có, phức tạp hơn. Thuật toán được đề xuất cho α-Neighbour p-Center là metaheuristic đầu tiên cho bài toán này.




## 6. Kết luận

Literature review này đã trình bày tổng quan về các phương pháp giải quyết bài toán P-Center Problem, tập trung vào các phương pháp SAT encoding, ILP/MIP, CP và các heuristic/metaheuristic. Mỗi phương pháp có những ưu điểm và nhược điểm riêng, phù hợp với các loại bài toán và kích thước dữ liệu khác nhau.

*   **SAT Encoding:** Phương pháp này, như được minh họa bởi Liu et al. [1], cho thấy tiềm năng trong việc tìm kiếm các giải pháp tối ưu hoặc gần tối ưu cho các bài toán P-Center bằng cách chuyển đổi chúng thành các bài toán con Set Covering và sử dụng các bộ giải SAT hiện đại. Ưu điểm chính là khả năng tận dụng sự phát triển nhanh chóng của các bộ giải SAT và khả năng giải quyết song song các bài toán con.

*   **ILP/MIP:** Các công thức lập trình nguyên tuyến tính hỗn hợp (MILP) như được đề xuất bởi Ales và Elloumi [2] cung cấp các mô hình chặt chẽ và có thể được giải quyết bằng các bộ giải tối ưu hóa thương mại mạnh mẽ như CPLEX. Các cải tiến trong công thức có thể giảm đáng kể số lượng ràng buộc và biến, dẫn đến thời gian giải quyết nhanh hơn, đặc biệt cho các trường hợp lớn hơn.

*   **CP (Constraint Programming):** Phương pháp lập trình ràng buộc, thường kết hợp với các heuristic và tìm kiếm cục bộ, như được trình bày bởi Iosif et al. [3], cho thấy sự linh hoạt trong việc xử lý các ràng buộc phức tạp và có thể hiệu quả cho các bài toán có nhiều giải pháp. Mặc dù các bộ giải CP không hoàn chỉnh có thể không đảm bảo tính tối ưu, chúng có thể tìm thấy các giải pháp chất lượng cao trong thời gian hợp lý cho các bài toán lớn.

*   **Heuristic và Metaheuristic:** Các phương pháp heuristic và metaheuristic, bao gồm Local Search, Tabu Search, Variable Neighborhood Search, Memetic Genetic Algorithm, và các biến thể khác như được Mousavi [4] nghiên cứu, là rất quan trọng cho các bài toán P-Center kích thước lớn, nơi các phương pháp chính xác trở nên không khả thi. Chúng cung cấp các giải pháp chất lượng tốt trong thời gian chấp nhận được, mặc dù không đảm bảo tính tối ưu.

Nhìn chung, việc lựa chọn phương pháp phụ thuộc vào yêu cầu cụ thể của bài toán, bao gồm kích thước, độ phức tạp của ràng buộc, và yêu cầu về tính tối ưu. Sự kết hợp giữa các kỹ thuật từ các lĩnh vực khác nhau (ví dụ: SAT với các quy tắc giảm dữ liệu, CP với tìm kiếm cục bộ) cho thấy một hướng đi đầy hứa hẹn để giải quyết hiệu quả hơn các bài toán P-Center Problem.

## 7. Tài liệu tham khảo

[1] X. Liu et al., “Effective Approaches to Solve P-Center Problem via Set Covering and SAT,” *IEEE Access*, vol. 8, pp. 161232–161245, 2020. [https://ieeexplore.ieee.org/document/9181511](https://ieeexplore.ieee.org/document/9181511)

[2] Z. Ales and S. Elloumi, “Compact MILP formulations for the $p$-center problem,” *arXiv preprint arXiv:2302.04591*, 2023. [https://arxiv.org/abs/2302.04591](https://arxiv.org/abs/2302.04591)

[3] P. Iosif, N. Ploskas, K. Stergiou, and D. C. Tsouros, “A CP/LS Heuristic Method for Maxmin and Minmax Location Problems with Distance Constraints,” in *30th International Conference on Principles and Practice of Constraint Programming (CP 2024)*, 2024, vol. 14, pp. 14:1–14:21. [https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CP.2024.14](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CP.2024.14)

[4] S. R. Mousavi, “Exploiting flat subspaces in local search for p-Center problem and two fault-tolerant variants,” *Computers & Operations Research*, vol. 149, p. 106023, 2023. [https://www.sciencedirect.com/science/article/pii/S0305054822002532](https://www.sciencedirect.com/science/article/pii/S0305054822002532)


