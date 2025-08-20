
import abcEconomics as abce
import random

class Firm(abce.Agent):
    """
    기업(Firm) 에이전트입니다.
    이 에이전트는 가계로부터 노동(labor)을 고용하여 상품(bread)을 생산하고,
    생산된 상품을 가계에 판매하여 이윤을 극대화하는 것을 목표로 합니다.
    """
    def init(self, num_households):
        """
        기업 에이전트의 초기 상태를 설정합니다.
        """
        # 거래할 가계의 수를 저장합니다.
        self.num_households = num_households

        # 생산된 빵(bread)을 저장할 공간을 만듭니다. 처음엔 0개입니다.
        self.create('bread', 0)
        
        # 기업 운영을 위한 초기 자본(money)을 보유합니다.
        self.create('money', 100)

        # 콥-더글라스 생산함수 Y = K^alpha * L^(1-alpha) 의 파라미터들
        # K (자본)는 단순화를 위해 1로 고정합니다.
        self.capital = 1
        self.alpha = 0.5  # 자본 소득 분배율

        # 지난 라운드의 수익을 기록하여 가격 결정에 활용합니다.
        self.last_revenue = 0
        
        # 이력(log)을 기록하여 시간 경과에 따른 변화를 관찰합니다.
        self.log_data = {'price': [], 'wage': [], 'production': [], 'profit': []}


    def hire_labor(self):
        """
        노동 시장에서 노동자를 고용합니다.
        랜덤한 가계 한 명을 선택하여 채용을 제안합니다.
        """
        # 고정된 임금(wage)으로 채용 공고(구매 오퍼)를 올립니다.
        wage = 10 
        quantity_to_hire = 1 # 1명의 노동자만 고용한다고 가정
        
        # 0부터 (총 가계 수 - 1) 사이에서 무작위 ID를 선택
        if self.num_households > 0:
            random_household_id = random.randint(0, self.num_households - 1)
            receiver = ('household', random_household_id)
            self.buy(receiver, good='labor', quantity=quantity_to_hire, price=wage)
        
        self.log_data['wage'].append(wage)


    def produce_bread(self):
        """
        고용한 노동(labor)을 사용하여 빵(bread)을 생산합니다.
        생산량은 콥-더글라스 생산함수에 의해 결정됩니다.
        """
        import math
        
        # 이번 라운드에 고용한 노동력의 양을 확인합니다.
        labor_hired = self.possession('labor')

        if labor_hired > 0:
            # 콥-더글라스 생산함수: Y = K^alpha * L^(1-alpha)
            production_qty = math.pow(self.capital, self.alpha) * math.pow(labor_hired, 1 - self.alpha)
            
            # 생산된 빵을 기업의 재고에 추가합니다.
            self.create('bread', production_qty)
            
            # 고용된 노동은 생산 과정에서 소모됩니다.
            self.destroy('labor', labor_hired)

            self.log_data['production'].append(production_qty)
            # print(f"{self.name} produced {production_qty:.2f} bread with {labor_hired} labor.")
        else:
            self.log_data['production'].append(0)


    def sell_bread(self):
        """
        생산한 빵을 상품 시장에 판매합니다.
        랜덤한 가계 한 명을 선택하여 판매를 제안합니다.
        """
        # 판매할 빵이 있다면
        if self.possession('bread') > 0:
            # 고정된 가격(price)으로 판매 오퍼를 올립니다.
            price = 12
            quantity_to_sell = self.possession('bread')
            
            if self.num_households > 0:
                random_household_id = random.randint(0, self.num_households - 1)
                receiver = ('household', random_household_id)
                self.sell(receiver, good='bread', quantity=quantity_to_sell, price=price)

            self.log_data['price'].append(price)


    def calculate_profit(self):
        """
        라운드가 끝난 후, 수익을 계산하고 이윤을 반환합니다.
        """
        # 현재 보유한 돈에서 지난 라운드 시작 시점의 돈을 빼서 수익을 계산합니다.
        # (더 정확하려면, 임금 지불 전후의 돈을 비교해야 합니다)
        # 여기서는 단순하게, (이번 라운드 판매수익) - (이번 라운드 임금지출) 로 계산합니다.
        
        # 이번 라운드에 번 돈 (판매 수익)
        revenue = self.possession('money') - 100 # 초기자본 제외
        
        # 이번 라운드에 쓴 돈 (임금)
        wage_paid = self.log_data['wage'][-1] * self.log_data['production'][-1] if self.log_data['production'] and self.log_data['wage'] else 0

        profit = revenue - wage_paid
        self.log_data['profit'].append(profit)

        # 다음 라운드를 위해 돈을 다시 초기 자본으로 리셋합니다.
        # (이 부분은 모델 설계에 따라 달라질 수 있습니다. 이윤을 누적시킬 수도 있습니다.)
        self.destroy('money', self.possession('money'))
        self.create('money', 100)
        
        return profit if profit > 0 else 0
