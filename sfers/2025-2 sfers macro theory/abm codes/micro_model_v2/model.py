
import abcEconomics as abce
from household import Household
from firm import Firm

class MarketModel(abce.Simulation):
    """
    미시경제 순환 모델의 시뮬레이션을 정의하고 실행하는 클래스입니다.
    """
    def __init__(self, num_households, num_firms, rounds):
        """
        시뮬레이션 환경을 초기화합니다.

        Args:
            num_households (int): 시뮬레이션에 참여할 가계 에이전트의 수.
            num_firms (int): 시뮬레이션에 참여할 기업 에이전트의 수.
            rounds (int): 시뮬레이션 총 라운드 수.
        """
        super().__init__(name='MicroModel')

        self.num_households = num_households
        self.num_firms = num_firms
        self.rounds = rounds
        
        # 에이전트들을 생성합니다.
        # 가계(Household) 그룹을 만듭니다.
        # num_firms 인자를 전달하여, 각 가계가 기업 수를 알게 합니다.
        self.households = self.build_agents(Household, 'household', self.num_households, num_firms=self.num_firms)
        
        # 기업(Firm) 그룹을 만듭니다.
        self.firms = self.build_agents(Firm, 'firm', self.num_firms, num_households=self.num_households)

        


    def run_round(self):
        """
        시뮬레이션의 한 라운드(예: 하루)를 실행합니다.
        각 단계는 경제 활동의 인과관계를 고려하여 순서대로 진행됩니다.
        """
        # --- 1. 노동 시장 (Labor Market) ---
        # 기업이 먼저 채용 공고(buy offer for labor)를 냅니다.
        for firm in self.firms:
            firm.hire_labor()
        # 가계는 그 공고를 보고 자신의 노동력을 판매(sell labor)합니다.
        for household in self.households:
            household.sell_labor()

        # --- 2. 생산 (Production) ---
        # 기업은 고용한 노동력을 사용하여 빵을 생산합니다.
        for firm in self.firms:
            firm.produce_bread()

        # --- 3. 상품 시장 (Goods Market) ---
        # 기업이 생산한 빵을 시장에 판매(sell offer for bread)합니다.
        for firm in self.firms:
            firm.sell_bread()
        # 가계는 소득을 사용하여 빵을 구매(buy bread)합니다.
        for household in self.households:
            household.buy_bread()

        # --- 4. 정산 및 이윤 분배 (Accounting & Profit Distribution) ---
        # 기업으로부터 이윤을 수집합니다.
        total_profit = sum([firm.calculate_profit() for firm in self.firms])

        # 총 이윤을 모든 가계에게 균등하게 배당금으로 분배합니다.
        if self.num_households > 0:
            dividend_per_household = total_profit / self.num_households
            if dividend_per_household > 0:
                for household in self.households:
                    household.create('money', dividend_per_household)

        # --- 5. 상태 기록 (Logging) ---
        # 각 에이전트의 현재 상태를 기록합니다.
        for household in self.households:
            household.log_status()

    def run(self):
        """
        설정된 라운드 수만큼 시뮬레이션을 실행합니다.
        """
        for r in range(self.rounds):
            print(f"--- Round {r+1}/{self.rounds} ---")
            self.run_round()

    def finish(self):
        """
        시뮬레이션이 끝난 후, 결과를 집계하고 출력합니다.
        """
        # 예를 들어, 모든 가계의 평균 자산 변화를 출력할 수 있습니다.
        import matplotlib.pyplot as plt
        
        # 모든 가계의 돈(money) 로그를 합산
        all_money_logs = [h.log_data['money'] for h in self.households]
        # 라운드별 평균 계산
        avg_money = [sum(round_data) / len(round_data) for round_data in zip(*all_money_logs)]

        # 모든 가계의 효용(utility) 로그를 합산
        all_utility_logs = [h.log_data['utility'] for h in self.households]
        avg_utility = [sum(round_data) / len(round_data) for round_data in zip(*all_utility_logs)]

        # 그래프 그리기
        plt.figure(figsize=(12, 5))

        plt.subplot(1, 2, 1)
        plt.plot(avg_money)
        plt.title("Average Household Money over Time")
        plt.xlabel("Round")
        plt.ylabel("Average Money")

        plt.subplot(1, 2, 2)
        plt.plot(avg_utility)
        plt.title("Average Household Utility over Time")
        plt.xlabel("Round")
        plt.ylabel("Average Utility")

        plt.tight_layout()
        plt.savefig("simulation_results.png")
        print("Simulation results saved to simulation_results.png")

