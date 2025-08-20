
# Micro-Foundation을 갖는 ABM 개발 튜토리얼 (with BeforeIT.jl)

이 문서는 Julia 언어와 `BeforeIT.jl` 패키지를 사용하여, 각 경제주체(agent)가 자신의 효용(utility)을 극대화하는 합리적인 의사결정을 내리는 간단한 Agent-Based Model (ABM)을 만드는 과정을 상세히 설명합니다.

**독자 대상**: Python 등 다른 언어에 경험이 있지만 Julia는 처음인 분들을 대상으로 합니다.

**핵심 목표**:
1.  Julia 프로젝트 환경 설정하기
2.  `BeforeIT.jl`을 이용해 ABM의 기본 구조 만들기
3.  `Optim.jl`을 이용해 Agent의 효용 극대화 문제 풀기
4.  시뮬레이션을 실행하고 결과 확인하기

---

## 1. 프로젝트 환경 설정

가장 먼저, 이 프로젝트에 필요한 패키지들을 설치하고 관리하기 위한 독립적인 환경을 설정합니다. Python의 `venv`와 유사한 기능입니다.

**단계 1: Julia 실행 및 폴더 이동**

먼저 Julia REPL(커맨드 창)을 실행하고, 아래 명령어를 입력하여 이 파일이 있는 폴더로 이동합니다.

```{julia}
cd("E:/github_jh/jh_hub/sfers/2025-2 sfers macro theory/abm codes/micro_foundation_tutorial")
```

**단계 2: 패키지 관리 모드 진입 및 환경 활성화**

Julia REPL에서 `]` 키를 누르면, 프롬프트가 `julia>`에서 `(@v1.10) pkg>` 와 같이 바뀌며 패키지 관리 모드로 들어갑니다. 여기서 다음 명령어를 입력하여 현재 폴더를 이 프로젝트의 독립적인 환경으로 지정합니다.

```pkg
activate .
```
이 명령은 현재 폴더에 `Project.toml`과 `Manifest.toml`이라는 파일을 생성하여, 이 프로젝트가 사용하는 패키지와 정확한 버전을 관리하게 해줍니다. Python의 `requirements.txt`와 유사하지만 훨씬 더 정교합니다.

**단계 3: 필수 패키지 추가**

이제 이 프로젝트에서 사용할 패키지들을 추가합니다.

1.  **`BeforeIT.jl`**: ABM의 기본 프레임워크입니다. 로컬에 있는 버전을 사용하도록 경로를 지정해줍니다.
2.  **`Optim.jl`**: 수치 최적화(numerical optimization)를 위한 패키지입니다. Agent의 효용 극대화 문제를 푸는 데 사용됩니다.
3.  **`Plots.jl`**: 시뮬레이션 결과를 시각화하기 위한 패키지입니다.
4.  **`Revise.jl`**: (강력 추천) 코드를 수정했을 때 Julia 세션을 재시작하지 않고도 변경사항을 즉시 반영해주는 필수 패키지입니다.

패키지 관리 모드(`pkg>`)에서 다음 명령어들을 차례로 실행하세요.

```pkg
add "E:\github_jh\jh_hub\sfers\2025-2 sfers macro theory\abm codes\BeforeIT.jl-main"
add Optim
add Plots
add Revise
```

이제 `Esc` 키나 `Backspace` 키를 눌러 다시 `julia>` 프롬프트로 돌아옵니다.

---

## 2. 경제 모형 설계

우리는 매우 간단한 경제를 가정합니다. 이 경제에는 **가계(Household)** 라는 단일 종류의 Agent만 존재합니다.

### 가계(Household)의 문제

가계는 **소비(consumption, c)** 와 **노동(labor, l)** 을 통해 효용을 얻습니다. 가계의 목표는 주어진 임금(wage, w) 하에서 자신의 효용을 극대화하는 것입니다.

-   **효용 함수 (Utility Function)**: 가계의 효용은 소비에서 오는 만족감과 노동에서 오는 비효용(disutility)으로 결정됩니다. 여기서는 로그-선형(log-linear) 형태의 효용 함수를 가정합니다.

    $
    \max_{c, l} \quad U(c, l) = \log(c) - \frac{l^2}{2}
    $

-   **예산 제약 (Budget Constraint)**: 가계의 소비는 노동 소득을 초과할 수 없습니다. 물가(price, p)가 1이라고 가정하면, 예산 제약은 다음과 같습니다.

    $
    c \le w \cdot l
    $

가계는 합리적이므로 항상 소득을 전부 소비합니다 (`c = w * l`). 따라서 가계의 문제는 노동(l)을 선택하여 효용을 극대화하는 단일 변수 문제로 단순화됩니다.

$
\max_{l} \quad \log(w \cdot l) - \frac{l^2}{2}
$

우리는 이 최적화 문제를 `Optim.jl`을 사용하여 수치적으로 풀 것입니다.

---

## 3. Julia 코드 구현

이제 실제 코드를 작성해 보겠습니다. 아래의 코드 블록들을 순서대로 Julia REPL에 복사하여 실행하거나, 하나의 `.jl` 파일로 만들어 실행할 수 있습니다.

### 3.1. 패키지 로딩 및 기본 설정

Python의 `import` 와 유사하게, Julia에서는 `using` 키워드로 패키지를 로드합니다.

```julia
# Revise.jl은 코드 수정 시 자동 리로딩을 위해 가장 먼저 로드하는 것이 좋습니다.
using Revise

# 필수 패키지들을 로드합니다.
using BeforeIT
using Optim
using Plots

# 플로팅 백엔드를 설정합니다. (처음 실행 시 약간의 시간이 걸릴 수 있습니다)
gr()
```

### 3.2. Agent 정의하기

`BeforeIT.jl`의 `Agent`를 상속받아 우리 모델의 `Household` Agent를 정의합니다.

-   `mutable struct`: Python의 `class`와 유사하게, 내장된 속성(field)들의 값을 변경할 수 있는 복합 데이터 타입입니다. (`struct`는 한 번 생성되면 값을 바꿀 수 없는 불변(immutable) 타입입니다.)
-   `<: Agent`: `Agent` 타입을 상속받는다는 의미입니다. `BeforeIT.jl`의 모든 Agent는 이를 따라야 합니다.
-   `id::Int`, `pos::Tuple{Int,Int}`: `BeforeIT.jl`이 요구하는 필수 속성입니다. `id`는 고유 번호, `pos`는 공간 모델을 위한 좌표(여기서는 사용하지 않음)입니다.
-   `labor::Float64`, `consumption::Float64`: 우리 모델에서 각 가계의 노동 공급량과 소비량을 저장하기 위한 속성입니다. `::Float64`는 이 속성이 64비트 부동소수점 숫자임을 명시하는 타입 선언입니다.

```julia
mutable struct Household <: Agent
    # BeforeIT.jl 필수 속성
    id::Int
    pos::Tuple{Int, Int}

    # 모델 고유 속성
    labor::Float64
    consumption::Float64
end
```

### 3.3. 모델(거시 환경) 정의하기

Agent들이 상호작용하는 거시적 환경인 `Model`을 정의합니다.

-   `<: Model`: `BeforeIT.jl`의 `Model` 타입을 상속받습니다.
-   `agents::Vector{Household}`: 모델에 포함된 모든 `Household` Agent들을 담는 벡터(리스트)입니다. `Vector{Household}`는 `Household` 타입의 객체만 담을 수 있음을 명시합니다.
-   `properties::Dict`: 모델의 거시 변수들을 저장하는 딕셔너리입니다. 여기서는 임금(`wage`)을 저장합니다. Python의 `dict`와 같습니다.

```julia
mutable struct EcoModel <: Model
    # BeforeIT.jl 필수 속성
    agents::Vector{Household}
    properties::Dict
end
```

### 3.4. Agent의 행동 함수 정의하기

이 부분이 모델의 핵심입니다. 각 `Household` Agent가 매 스텝마다 어떻게 행동할지를 정의하는 함수입니다.

-   `agent_step!(h::Household, m::EcoModel)`: Julia의 가장 강력한 기능인 **다중 디스패치(Multiple Dispatch)** 를 보여줍니다. 이 함수는 첫 번째 인자가 `Household` 타입이고 두 번째 인자가 `EcoModel` 타입일 때만 호출됩니다. 함수 이름에 `!`를 붙이는 것은 함수가 인자의 내용을 변경할 수 있다는 Julia의 관례입니다.

```julia
function agent_step!(h::Household, m::EcoModel)
    # 1. 모델의 거시 변수(임금)를 가져옵니다.
    wage = m.properties[:wage]

    # 2. 최적화할 목적 함수를 정의합니다.
    #    l은 벡터(l[1])로 전달됩니다. Optim.jl의 기본 형식입니다.
    #    Optim.jl은 기본적으로 최소화를 하므로, 효용 극대화를 위해 효용 함수에 -1을 곱해줍니다.
    objective_function(l) = -(log(wage * l[1]) - (l[1]^2) / 2)

    # 3. 노동 공급량의 초기 추정치를 설정합니다.
    initial_l = [1.0]

    # 4. 제약 없는 최적화 문제를 풉니다.
    #    Optim.optimize(목적함수, 초기값, 최적화알고리즘)
    result = Optim.optimize(objective_function, initial_l, NelderMead())

    # 5. 최적화 결과를 Agent의 속성에 저장합니다.
    #    Optim.minimizer(result)가 최적해(노동량)를 반환합니다.
    optimal_labor = Optim.minimizer(result)[1]
    h.labor = optimal_labor
    h.consumption = wage * optimal_labor
end
```

### 3.5. 시뮬레이션 초기화 및 실행

이제 모델을 만들고, Agent를 추가하고, 시뮬레이션을 실행할 차례입니다.

```julia
# 1. 모델 초기화
#    - 100명의 가계(Household)를 생성합니다.
#    - 초기 임금(wage)은 1.5로 설정합니다.
num_agents = 100
model = EcoModel(
    [Household(i, (0,0), 0.0, 0.0) for i in 1:num_agents],
    Dict(:wage => 1.5)
)

# 2. 데이터 수집 설정
#    - `adata`는 각 스텝마다 수집할 agent 데이터의 리스트입니다.
#    - `[:labor, :consumption]`은 각 agent의 labor와 consumption 속성을 수집하라는 의미입니다.
#    - `agg`는 수집한 데이터를 어떻게 집계할지를 정의합니다. 여기서는 평균(mean)을 계산합니다.
adata = [(:labor, agg=mean), (:consumption, agg=mean)]

# 3. 시뮬레이션 실행
#    - run!(모델, agent_step_함수, 총_스텝_수, adata=수집할_agent_데이터)
#    - 이 시뮬레이션은 1 스텝만 실행합니다. (임금이 고정되어 있으므로 여러 스텝을 돌려도 결과는 같습니다)
results, _ = run!(model, agent_step!, 1, adata=adata)

# 4. 결과 출력
println(results)
```

### 3.6. (심화) 임금이 변하는 시나리오

이번에는 매 스텝마다 임금이 조금씩 오르는 동적인 시나리오를 시뮬레이션해 보겠습니다. 이를 위해 `model_step!` 함수를 정의합니다. `model_step!`은 모든 agent들의 `agent_step!`이 끝난 후 각 스텝마다 한 번씩 호출됩니다.

```julia
# 모델 스텝 함수: 매 스텝마다 임금을 0.1씩 올립니다.
function model_step!(m::EcoModel)
    m.properties[:wage] += 0.1
end

# 모델 재초기화
model = EcoModel(
    [Household(i, (0,0), 0.0, 0.0) for i in 1:num_agents],
    Dict(:wage => 1.0) # 시작 임금을 1.0으로 설정
)

# 20 스텝 동안 시뮬레이션 실행 (agent_step!과 model_step! 모두 사용)
results_dynamic, _ = run!(model, agent_step!, model_step!, 20, adata=adata)

# 결과 출력
println(results_dynamic)

# 결과 시각화
plot(results_dynamic.mean_labor, label="Average Labor Supply", xlabel="Time Step", ylabel="Labor")
plot!(twinx(), results_dynamic.mean_consumption, label="Average Consumption", ylabel="Consumption", color=:red, legend=:bottomright)

```

위 코드를 실행하면 임금이 상승함에 따라 가계들이 노동 공급과 소비를 늘리는 것을 보여주는 그래프가 나타날 것입니다. 이것이 바로 미시적 최적화 원리가 거시적 현상으로 나타나는 ABM의 특징입니다.

---

## 4. 결론 및 다음 단계

이 튜토리얼을 통해 Julia와 `BeforeIT.jl`, `Optim.jl`을 사용하여 micro-foundation을 갖는 간단한 ABM을 구축하는 전 과정을 살펴보았습니다.

**다음 단계 아이디어**:
-   **기업(Firm) Agent 추가**: 이윤을 극대화하는 기업 Agent를 추가하여 노동 수요를 내생적으로 결정하게 할 수 있습니다.
-   **시장 메커니즘 구현**: `model_step!`에서 노동 시장과 상품 시장의 수요-공급을 일치시켜 임금과 물가가 동적으로 결정되도록 모델을 확장할 수 있습니다.
-   **다양한 효용/생산 함수 실험**: 더 복잡하고 현실적인 함수들을 사용하여 모델을 정교화할 수 있습니다.

이 파일의 코드들을 Julia REPL에 단계적으로 실행해보면서 각 부분이 어떻게 작동하는지 직접 확인해 보시기를 권장합니다.
